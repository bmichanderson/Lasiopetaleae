#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: June 2023
# Description: stitch selected blast hit contigs together to make supercontigs
# Note: it requires each sample to have 1) a contigs file (named sample.fa)
#	and 2) a file of the selected blast hits (sample_blast_out.tab) with specific headers:
#	sacc slen pident length qstart qend sstart send
#	(in the working directory where this script is called)
##########################


import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to stitch ' +
	'selected blast hit contigs together, separated by 100 Ns, into supercontigs')


# add arguments to parse
parser.add_argument('con_fastas', type=str, help='The contig fasta files', nargs='*')
parser.add_argument('-s', type = str, dest = 'samples_file',
	help = 'A file with sample names, one per line')
parser.add_argument('-v', type = str, dest = 'overlap',
	help = 'Allowed overlap of hits (bp) to keep them [default 20]; otherwise, the best hit (pident * length) is kept')

# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
con_fastas = args.con_fastas
samples_file = args.samples_file
overlap = args.overlap

if not samples_file:
	print('Please specify a samples file with -s \n')
	parser.print_help(sys.stderr)
	sys.exit(1)

if not con_fastas:
	print('Please specify contig fasta files (*.fa)\n')
	parser.print_help(sys.stderr)
	sys.exit(1)

if not overlap:
	overlap = 20

fasta_files = []
for file in con_fastas:
	fasta_files.append(file)


# read in the samples file and extract names
samples = []
with open(samples_file, 'r') as infile:
	for line in infile:
		samples.append(str(line.strip()))


# cycle through the samples
# for each, evaluate the blast hits and note contigs
# then go through the contigs file and extract them
# pasting them together into an output file
master_drop_count = []
for sample in samples:
	blast_file = str(sample) + '_blast_out.tab'
	contig_file = [i for i in fasta_files if (str(sample) + '.fa') in i][0]

	# open the blast file and iterate through the hits,
	# storing the contigs and order
	blastdf = pd.read_csv(blast_file, sep = '\t', header = 0)
	sorted_df = blastdf.sort_values('qstart', ascending = True)
	sorted_df.reset_index(drop = True, inplace = True)
	contigs = []
	orientation = []
	drop_contigs = []
	if len(set(sorted_df['sacc'])) == 1:	# a single contig (may have multiple hits)
		contigs.append(sorted_df['sacc'][0])
		if sorted_df['sstart'][0] > sorted_df['send'][0]:
			this_orient = 'minus'
		else:
			this_orient = 'plus'
		orientation.append(this_orient)
	else:	# multiple hits, which are ordered from the start of the query
			# if hits do not overlap, keep the first contig and continue
			# if hits overlap, keep the greater of pid*hitlen
		max_coord = 0
		for row in range(len(sorted_df) - 1):
			this_contig = sorted_df['sacc'][row]
			this_pidxlen = sorted_df['pident'][row] * sorted_df['length'][row]
			this_coord = (sorted_df['qstart'][row], sorted_df['qend'][row])
			if sorted_df['sstart'][row] > sorted_df['send'][row]:
				this_orient = 'minus'
			else:
				this_orient = 'plus'

			next_contig = sorted_df['sacc'][row + 1]
			next_pidxlen = sorted_df['pident'][row + 1] * sorted_df['length'][row + 1]
			next_coord = (sorted_df['qstart'][row + 1], sorted_df['qend'][row + 1])
			if sorted_df['sstart'][row + 1] > sorted_df['send'][row + 1]:
				next_orient = 'minus'
			else:
				next_orient = 'plus'

			if this_coord[1] > max_coord:			# to avoid internal hits
				if all([next_contig in contigs, this_contig not in contigs]):	# internal to a contig already selected
					drop_contigs.append(this_contig)

					continue
				else:
					max_coord = this_coord[1]
			else:	# internal hit of the previous contig
				if this_contig not in contigs:
					drop_contigs.append(this_contig)

				continue

			if this_coord[1] < next_coord[0] + int(overlap):		# does not overlap by more than specified amount
				if all([this_contig not in contigs, this_contig not in drop_contigs]):
					contigs.append(this_contig)
					orientation.append(this_orient)

				if all([row == len(sorted_df) - 2, next_contig != this_contig,
						next_contig not in contigs, next_contig not in drop_contigs]):
					contigs.append(next_contig)
					orientation.append(next_orient)

			else:										# overlaps
				if this_pidxlen >= next_pidxlen:		# if the first hit is "better"
					if all([this_contig not in contigs, this_contig not in drop_contigs]):
						contigs.append(this_contig)
						orientation.append(this_orient)
					elif all([this_contig in drop_contigs, row == len(sorted_df) - 2,
							next_contig != this_contig, next_contig not in contigs,
							next_contig not in drop_contigs]):
						contigs.append(next_contig)
						orientation.append(next_orient)

					if all([next_contig != this_contig, next_contig not in contigs]):
						drop_contigs.append(next_contig)

				else:									# if the next hit is "better"
					if all([next_contig not in contigs, next_contig not in drop_contigs]):
						contigs.append(next_contig)
						orientation.append(next_orient)
					elif all([next_contig in drop_contigs, this_contig not in contigs,
							this_contig not in drop_contigs]):
						contigs.append(this_contig)
						orientation.append(this_orient)

					if all([next_contig != this_contig, this_contig not in contigs]):
						drop_contigs.append(this_contig)


	if len(list(set(contigs) & set(drop_contigs))) > 0:		# if there is a contig supposed to be dropped but kept
		print('A contig is being kept for sample ' + sample + ' that should not!')


	# now create the output for the sample, concatenating sequences
	# where required, separated by 100 Ns
	contig_fastas = []
	with open(contig_file, 'r') as infile:
		entries = SeqIO.parse(infile, 'fasta')
		for entry in entries:
			contig_fastas.append(entry)

	with open(str(sample) + '_stitched.fasta', 'w') as outfile:
		myid = str(sample)
		mydescription = str(sample)
		seqs = []
		for index, contig in enumerate(contigs):
			for entry in contig_fastas:
				if entry.id == contig:
					if orientation[index] == 'minus':
						this_seq = entry.seq.reverse_complement()
					else:
						this_seq = entry.seq
					seqs.append(this_seq)
					break
		myseq = Seq((100 * 'N').join([str(seq) for seq in seqs]))
		myrec = SeqRecord(myseq, id = myid, description = mydescription)
		SeqIO.write(myrec, outfile, 'fasta')

	# if there were dropped contigs, add them to the tally
	if len(set(drop_contigs)) > 0:
		master_drop_count.append([sample, list(set(drop_contigs))])

# report dropped contigs, if present
if len(master_drop_count) > 0:
	print('\nAdditional hits present:')
	for entry in master_drop_count:
		sample = str(entry[0])
		drop_list = entry[1]
		print('\n' + sample + '\t' + str(len(drop_list)) + '\t' +
			', '.join([str(i) for i in drop_list]) + '\n')
