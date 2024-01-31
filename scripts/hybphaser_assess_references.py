#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: April 2023
# Description: using the results of clade association, choose the references to phase to
# Note: Arguments are the clade association results, a csv file of clade references used for clade association,
#		the proportion of the top mapping association value for a reference to be kept (default: 0.33),
#       and an optional file of sample numbers of hybrids that need to be phased, one per line
##########################


import sys
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to assess the results of clade association in HybPhaser')


# add arguments to parse
parser.add_argument('-c', type = str, dest = 'clade_assoc', help = 'The clade association results file from HybPhaser')
parser.add_argument('-r', type = str, dest = 'clade_refs', help = 'The csv file of clade references (including abbreviations)')
parser.add_argument('-s', type = str, dest = 'hybrids', help = 'An optional file with sample numbers to be phased (one per line)')
parser.add_argument('-p', type = str, dest = 'propor', help = 'The minimum proportion of the top value to keep a reference (default: 0.33)')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
clade_assoc = args.clade_assoc
clade_refs = args.clade_refs
hybrids = args.hybrids
propor = args.propor

if not clade_assoc:
	print('Please provide the HybPhaser clade association results\n')
	parser.print_help(sys.stderr)
	sys.exit(1)
elif not clade_refs:
	print('Please provide the csv file of clade references used by HybPhaser\n')
	parser.print_help(sys.stderr)
	sys.exit(1)

if not hybrids:
	print('No subset of samples to evaluate specified, so evaluating all samples\n')

if not propor:
	propor = 0.33
else:
	propor = float(propor)


# read in the results of clade association and capture relevant info
assoc_list = []
samples = []
with open(clade_assoc, 'r') as infile:
	for lineno, line in enumerate(infile):
		pieces = line.strip().split(',')
		if lineno == 0: 	# the first line is a header with refrence abbreviations
			ref_list = [item.strip('"') for item in pieces[1: ]]
		else:
			sampleid = str(pieces[0].strip('"'))
			val_list = [float(item) for item in pieces[1: ]]
			assoc_list.append((sampleid, val_list))
			samples.append(sampleid)

print('Read in results for ' + str(len(samples)) + 
	' samples associated with ' + str(len(ref_list)) + ' references')


# read in the clade references and make a dictionary
# the references are in the form: sample,abbreviation
# we want to use the abbreviation as the key, as that's what's used in the results file
refdict = {}
with open(clade_refs, 'r') as infile:
	for lineno, line in enumerate(infile):
		if lineno == 0:		# the first line is a header
			continue
		else:
			pieces = line.strip().split(',')
			refdict[str(pieces[1])] = str(pieces[0])


# read in the hybrids list, if present
assess_samples = []
if hybrids:
	with open(hybrids, 'r') as infile:
		for line in infile:
			assess_samples.append(str(line.strip()))
	
	print('Read in ' + str(len(assess_samples)) + ' samples to assess')
else:
	assess_samples = samples


# iterate through the samples, choosing references
# output the samples to phase and the references in the csv format:
# sample,ref1,ref1abb,ref2,ref2abb,ref3,ref3abb,ref4,ref4abb
with open('phasing_prep.csv' , 'w') as outfile:
	outfile.write('sample,ref1,ref1abb,ref2,ref2abb,ref3,ref3abb,ref4,ref4abb\n')
	count_samples = 0
	for sample in assess_samples:
		found = False
		needed = True
		refs = []
		outlist = []
		for entry in assoc_list:
			if entry[0] == sample:
				valsum = sum(entry[1])		# sum the associations for all refs
				# sort the values in descending order; grab the indices
				indices = sorted(range(len(entry[1])), key = entry[1].__getitem__, reverse = True)
				maxindex = indices[0]
				maxval = entry[1][maxindex]
				if (maxval / valsum) > 0.80:	# if the max association is more than 80%
					needed = False				# no need to phase
					break
				else:
					# keep the top two at least
					refs.append(ref_list[indices[0]])
					refs.append(ref_list[indices[1]])
					for index in indices[2: ]:
						if (entry[1][index] / maxval) >= propor:		# if an association is at least x of top
							refs.append(ref_list[index])				# keep it
						else:
							break

				found = True
				break

		if not found:
			if not needed:
				print('Sample ' + str(sample) + ' not requiring phasing (> 80% to one ref)')
			else:
				print('\nSample ' + str(sample) + ' not found in clade association results\n')
		else:
			for ref in refs[0: 4]:		# output the top 4
				outlist.append(refdict[ref] + ',' + ref)
			outfile.write(sample + ',' + ','.join(outlist) + '\n')
			count_samples = count_samples + 1

# report completion
print('Wrote selected references for ' + str(count_samples) + ' samples to \"phasing_prep.csv\"')
