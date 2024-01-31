#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: April 2023
# Modified: June 2023
# Description: create loci multifastas from input sample multifastas with entries ">taxon-locus"  
# Note: it will create a file for each locus (named with the locus name) in the current directory
#	It assumes each input fasta file is for a single sample (not mixed)
##########################


import sys
import argparse
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to create loci multifastas from sample multifastas')


# add arguments to parse
parser.add_argument('fastas', type=str, help='The sample fasta files', nargs='*')
parser.add_argument('-s', type = str, dest = 'samples_file', help = 'A file with sampleIDs, one per line')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
fastas = args.fastas
samples_file = args.samples_file

if not samples_file:
	print('Samples file not specified, so will include all samples\n')

if not fastas:
	print('Please supply at least one sample fasta file\n')
	parser.print_help(sys.stderr)
	sys.exit(1)


# if there is a samples file, read it
samples = []
if samples_file:
	with open(samples_file, 'r') as infile:
		for line in infile:
			samples.append(str(line.strip()))


# iterate through the sample fasta files and collect loci
locus_list = []
for fasta_file in fastas:
	with open(fasta_file, 'r') as infile:
		for fasta_entry in SeqIO.parse(infile, 'fasta'):
			pieces = fasta_entry.id.split('-')
			sample = pieces[0]
			locus = pieces[1]
			if not samples_file:
				if sample not in samples:
					samples.append(sample)
				else:
					print('Found a duplicate (sample ' + str(sample) + ') input file! Skipping!')
					break

				found = False
				for locus_entry in locus_list:
					if locus_entry[0] == locus:
						fasta_entry.id = sample
						fasta_entry.description = sample
						locus_entry[1].append(fasta_entry)
						found = True
						break

				if not found:
					fasta_entry.id = sample
					fasta_entry.description = sample
					locus_list.append((locus, [fasta_entry]))

			else:		# a samples list was provided
				if sample in samples:
					found = False
					for locus_entry in locus_list:
						if locus_entry[0] == locus:
							fasta_entry.id = sample
							fasta_entry.description = sample
							locus_entry[1].append(fasta_entry)
							found = True
							break

					if not found:
						fasta_entry.id = sample
						fasta_entry.description = sample
						locus_list.append((locus, [fasta_entry]))
				else:
					print('Skipping sample ' + str(sample) + ' as it was not in sample list')
					break


# output the loci files in the current directory, named by locus  
for locus_entry in locus_list:
	locus = locus_entry[0]
	fastas = locus_entry[1]
	with open(str(locus) + '.fasta', 'w') as outfile:
		SeqIO.write(fastas, outfile, 'fasta')

	print('Wrote ' + str(len(fastas)) + ' entries for locus ' + str(locus))
