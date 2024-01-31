#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: Nov 2022
# Modified: Mar 2023, Jun 2023
# Description: stitch selected blast hit contigs together to make supercontigs
# 	for each locus in a SECAPR run
# Note: it requires each sample (a number) to have 1) a contigs file (sample.fa)
#	and 2) a file of the selected blast hits from SECAPR find_target_contigs
#	(in the working directory where this script is called)
##########################


import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to stitch SECAPR ' +
	'selected blast hit contigs together, separated by 100 Ns, into supercontigs')


# add arguments to parse
parser.add_argument('con_fastas', type=str, help='The contig fasta files', nargs='*')
parser.add_argument('-s', type = str, dest = 'samples_file',
	help = 'A file with sample names, one per line')
parser.add_argument('-t', type = str, dest = 'trans_file',
	help = 'The locus name translation file (reference_fasta_header_info.txt)')
parser.add_argument('-v', type = str, dest = 'overlap',
	help = 'Allowed overlap of hits (bp) to keep them [default 20]; otherwise, the best hit (pident * length) is kept')

# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
con_fastas = args.con_fastas
samples_file = args.samples_file
trans_file = args.trans_file
overlap = args.overlap

if not samples_file:
	print('Please specify a samples file with -s \n')
	parser.print_help(sys.stderr)
	sys.exit(1)
elif not con_fastas:
	print('Please specify contig fasta files (*.fa)\n')
	parser.print_help(sys.stderr)
	sys.exit(1)

if not trans_file:
	print('No reference_fasta_header_info.txt supplied, so renaming by SECAPR default\n')

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


# read in the translation file (if present) and create dictionary
if trans_file:
	trans_dict = {}
	with open(trans_file, 'r') as infile:
		for line in infile:
			pieces = line.strip().split()
			entry = { str(pieces[0]): str(pieces[1]).split('-')[1]}
			trans_dict.update(entry)


# cycle through the samples
# for each, evaluate the blast hits and note contigs
# then go through the contigs file and extract them
# pasting them together into an output file
master_drop_count = []
for sample in samples:
	blast_file = str(sample) + '_selected_blast_hits.txt'
	contig_file = [i for i in fasta_files if (str(sample) + '.fa') in i][0]

	# open the blast file and iterate through the hits,
	# storing the contigs and order for each locus
	blastdf = pd.read_csv(blast_file, sep = '\t', header = 0)
	sorted_df = blastdf.sort_values(['qseqid', 'qstart'], ascending = [True, True])
	drop_count = []
	loci = []
	for locus in set(blastdf['qseqid']):
		contigs = []
		orientation = []
		drop_contigs = []
		mini_df = sorted_df[sorted_df['qseqid'] == locus]	# subset for that locus
		mini_df.reset_index(drop = True, inplace = True)
		if len(mini_df) == 1:	# a single hit
			contigs.append(mini_df['sseqid'][0])
			orientation.append(mini_df['sstrand'][0])
		else:	# multiple hits, which are ordered from the start of the query
				# if hits do not overlap, keep the first contig and continue
				# if hits overlap, keep the greater of pid*hitlen
			max_coord = 0
			for row in range(len(mini_df) - 1):
				this_contig = mini_df['sseqid'][row]
				this_pidxlen = mini_df['pident'][row] * mini_df['length'][row]
				this_coord = (mini_df['qstart'][row], mini_df['qend'][row])
				next_contig = mini_df['sseqid'][row + 1]
				next_pidxlen = mini_df['pident'][row + 1] * mini_df['length'][row + 1]
				next_coord = (mini_df['qstart'][row + 1], mini_df['qend'][row + 1])

				if this_coord[1] > max_coord:			# to avoid hits that don't overlap end
					max_coord = this_coord[1]
				else:
					if this_contig not in contigs:
						drop_contigs.append(this_contig)

					continue

				if this_coord[1] < next_coord[0] + int(overlap):		# does not overlap by more than specified amount
					if all([this_contig not in contigs, this_contig not in drop_contigs]):
						contigs.append(this_contig)
						orientation.append(mini_df['sstrand'][row])

					if all([row == len(mini_df) - 2, next_contig != this_contig, next_contig not in drop_contigs]):
						contigs.append(next_contig)
						orientation.append(mini_df['sstrand'][row + 1])

				else:										# overlaps
					if this_pidxlen >= next_pidxlen:		# if the first hit is "better"
						if all([this_contig not in contigs, this_contig not in drop_contigs]):
							contigs.append(this_contig)
							orientation.append(mini_df['sstrand'][row])
						elif all([this_contig in drop_contigs, row == len(mini_df) - 2,
								next_contig != this_contig, next_contig not in drop_contigs]):
							contigs.append(next_contig)
							orientation.append(mini_df['sstrand'][row + 1])

						if all([next_contig != this_contig, next_contig not in contigs]):
							drop_contigs.append(next_contig)

					else:									# if the next hit is "better"
						if all([next_contig not in contigs, next_contig not in drop_contigs]):
							contigs.append(next_contig)
							orientation.append(mini_df['sstrand'][row + 1])
						elif all([next_contig in drop_contigs, this_contig not in contigs, this_contig not in drop_contigs]):
							contigs.append(this_contig)
							orientation.append(mini_df['sstrand'][row])
						
						if all([next_contig != this_contig, this_contig not in contigs]):
							drop_contigs.append(this_contig)

		if len(list(set(contigs) & set(drop_contigs))) > 0:		# if there is a contig supposed to be dropped but kept
			print('A contig is being kept for sample ' + sample + ' and locus ' + str(locus) + ' that should not!')

		# add the list of contigs to the locus list
		loci.append([locus, contigs, orientation])
		# if there were droppped contigs, add them to the tally
		if len(set(drop_contigs)) > 0:
			drop_count.append([locus, list(set(drop_contigs))])

	# now create the output for the sample, concatenating sequences
	# where required, separated by 100 Ns
	contig_fastas = []
	with open(contig_file, 'r') as infile:
		entries = SeqIO.parse(infile, 'fasta')
		for entry in entries:
			contig_fastas.append(entry)

	with open(str(sample) + '_targetcons.fasta', 'w') as outfile:
		for locus_list in loci:
			locus = locus_list[0]
			contigs = locus_list[1]
			orientation = locus_list[2]
			if trans_file:
				myid = str(sample) + '-' + trans_dict[str(locus)]
				mydescription = myid
			else:
				myid = str(locus) + '_' + str(sample)
				mydescription = '|' + str(locus)
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

	# if there were dropped contigs, add them to the list
	if len(drop_count) > 0:
		master_drop_count.append([sample, drop_count])

# report dropped contigs, if present
if len(master_drop_count) > 0:
	print('\nPotential paralogs present:')
	with open('dropped_contigs.txt', 'w') as outfile:
		outfile.write('Sample\tLocus\tContigs\n')
		for entry in master_drop_count:
			print('\n' + str(entry[0]) + '_(' + str(len(entry[1])) + '):\t' +
				' '.join([str(i) for i in [x[0] for x in entry[1]]]))
			for subentry in entry[1]:
				outfile.write(str(entry[0]) + '\t' + str(subentry[0]) +
					'\t' + ','.join(subentry[1]) + '\n')
