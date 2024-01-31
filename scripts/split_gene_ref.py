#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: Nov 2022
# Modified: Feb 2023 (change output to keep gene info)
# Description: split genes from a single multifasta (arg) into individual files
# Note: it will create a file for each gene (named with the gene name) in the current directory
#		It expects the reference fasta file to have entries in the form: ">taxon-gene "
##########################


import sys
import argparse
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to split genes from a multifasta as individual files')


# add arguments to parse
parser.add_argument('multifasta', type=str, help='The multifasta')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
multifasta = args.multifasta


# read in the multifasta and parse the entry names
# writing new files named after the genes
counter = 0
with open(multifasta, 'r') as infile:
	entries = SeqIO.parse(infile, 'fasta')
	for entry in entries:
		gene = entry.id.split('-')[1]
		with open(str(gene) + '.fasta', 'w') as outfile:
			SeqIO.write(entry, outfile, 'fasta-2line')
			counter = counter + 1

# record completion
print('Split ' + str(counter) + ' genes from the reference into new files')
