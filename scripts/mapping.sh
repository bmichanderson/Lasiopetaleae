#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: Oct 2022
# Modified: Jun 2023
# Description: run mapping of Illumina paired-end reads on Gadi
# NOTE: this script could be called from a separate job using nci-parallel
#	Arguments are two read files and a reference fasta file
#	Read pair files should be named with the same first number before an underscore
######################


# set a PBS directive in case this is run independently (so args can be specified)
#PBS -v read1,read2,ref


# set max memory to specify (limit to 85% when passing to bbmap scripts)
# based on my set up, it will be 2 GB per cpu
if [ -z $cores_per ]; then
	# this was launched independently from the launch script
	cores_per="$PBS_NCPUS"
fi
max_mem=$(echo "${cores_per}*2*85/100" | bc)g
echo "This job is going to use $cores_per cores with max mem of $max_mem"


# define where bbmap is
bbmap_loc="/g/data/nm31/bin/bbmap"


# bbmap args
minid="0.95"		# minimum identity for read mapping
ambig="random"		# the way to treat ambiguous mapping; random = randomly assign read; best = first best location
pairlen="400"	# maximum allowed distance between paired reads; insert = pairlen + read1 + read2
pairedonly="t"	# whether to map paired reads only


# read and check args
if [ -n "$read1" ] && [ -n "$read2" ] && [ -n "$ref" ]; then
	read1="$(readlink -f $read1)"
	read2="$(readlink -f $read2)"
	ref="$(readlink -f $ref)"
elif [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo -e "\nPlease specify input reads and ref as arguments 1, 2 and 3\n"
	exit 1
else
	read1="$(readlink -f $1)"
	read2="$(readlink -f $2)"
	ref="$(readlink -f $3)"
fi


# change to the job directory then make a directory for results
# also start a log file based on the read file name
cd $PBS_JOBFS
file_name=$(basename ${read1})
prefix=${file_name/_*/}
if [ ! -d "$prefix" ]; then
	mkdir -p "$prefix"
fi
cd "$prefix"
logfile="$prefix".log


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting Illumina read mapping at $(date)\n**********" >> $logfile 2>&1


# load modules
module load java/jdk-8.40 >> $logfile 2>&1


# run read mapping
"$bbmap_loc"/bbmap.sh ref="$ref" in1="$read1" in2="$read2" outm=map.sam bs=bs.sh nodisk minid="$minid" \
	ambiguous="$ambig" pairedonly="$pairedonly" pairlen="$pairlen" threads="$cores_per" -Xmx"$max_mem" \
	>> $logfile 2>&1
./bs.sh >> $logfile 2>&1
rm bs.sh map.sam


# move the results back to where the job was launched from  
cd ..
rsync -rut "$prefix" $PBS_O_WORKDIR/


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished Illumina read mapping at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
echo "This job finished"
