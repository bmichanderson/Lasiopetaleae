#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: Jan 2023
# Description: format phased read files (*.fastq)  
# NOTE: submit this specifying the reads directory (default: current directory):
#	qsub -v reads_dir="____" hybphaser_phasedreadformat.sh
######################


# Set PBS directives

#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v reads_dir


# set the container (with BBMap in it)
container="/g/data/nm31/ben_anderson/singularity/hybphaser.sif"


# set up a master log file
logfile="$(pwd)"/hybphaser_phasedreadformat_"${PBS_JOBID/\.gadi-pbs/}".log


# Load the singularity module
module load singularity >> $logfile 2>&1


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nFormatting reads files at $(date)\n**********" >> $logfile 2>&1


# Check the argument input and change to the reads directory
if [ -z "$reads_dir" ]; then
	reads_dir="$(pwd)"
else
	reads_dir="$(readlink -f $reads_dir)"
fi

cd $reads_dir


# compress the reads with pigz
for file in *.fastq; do
	pigz "$file"
done


# rename the files
for file in *.fastq.gz; do
	mv "$file" "${file/_to_/}"
done


# run reformat.sh to change the interleaved reads into paired files  
# remove the interleaved file after  
for file in *.fastq.gz; do
	singularity exec -H "$(pwd)" "$container" reformat.sh in="$file" out1="${file/.fastq/_R1.fastq}" \
		out2="${file/.fastq/_R2.fastq}"  >> $logfile 2>&1
	if [ -f ${file/.fastq/_R1.fastq} ]; then
		rm $file
	fi
done


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished formatting reads files after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
