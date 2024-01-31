#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: April 2023
# Updated: May 2023 (for paralogs)
# Description: collect desired loci from HybPiper results files
# Note: Arguments are a file with tab-separated sampleID and locus name, one per line,
#		and the results fasta files (format: "locusname.fasta") from a HybPiper run
#		The fasta entries should be named as ">sampleID ..."
#		Alternatively, if these are paralogs files, entries can be named ">sampleID.main" or ">sampleID.0" etc.
##########################


import sys
import os
import argparse
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to collect desired loci from HybPiper results')


# add arguments to parse
parser.add_argument('fastas', type=str, help='The results loci fasta files', nargs='*')
parser.add_argument('-l', type = str, dest = 'list_file', help = 'A file with tab-separated sampleID and locus name, one per line')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
fastas = args.fastas
list_file = args.list_file

if not list_file:
	print('Please specify a file sampleID and locus names using \"-l\"\n')
	parser.print_help(sys.stderr)
	sys.exit(1)
elif not fastas:
	print('Please supply at least one result locus fasta file\n')
	parser.print_help(sys.stderr)
	sys.exit(1)


# read in the list file and capture into a list of loci and lists of samples
master_list = []
loci = []
with open(list_file, 'r') as infile:
	for line in infile:
		pieces = line.strip().split()
		sample = str(pieces[0])
		locus = str(pieces[1])
		if locus not in loci:
			loci.append(locus)
			master_list.append((locus, [sample]))
		else:
			for entry in master_list:
				if entry[0] == locus:
					entry[1].append(sample)
					break


# iterate through the input fasta files and output the desired entries
for fasta_file in fastas:
	filename = os.path.basename(fasta_file)
	locus = str(filename.split('.')[0])
	if locus not in loci:
		continue
	else:
		for entry in master_list:
			if entry[0] == locus:
				sample_list = entry[1]
				break
	
	with open(fasta_file, 'r') as infile:
		fasta_list = []
		for fasta_entry in SeqIO.parse(infile, 'fasta'):
			if fasta_entry.id in sample_list:
				fasta_entry.description = fasta_entry.id
				fasta_list.append(fasta_entry)
			elif fasta_entry.id.split('.')[0] in sample_list:		# in case these are paralogs
				fasta_entry.description = fasta_entry.id
				fasta_list.append(fasta_entry)
			else:
				continue

	with open('collected_' + filename, 'w') as outfile:
		SeqIO.write(fasta_list, outfile, 'fasta')
