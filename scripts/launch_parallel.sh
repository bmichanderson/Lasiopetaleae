#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: June 2023
# Description: submit a list of multiple jobs to be parallelised on Gadi
# Note: The `calls.txt` file should have a call to a job per line
#	The arguments for this script include the path to the calls file,
#	cores to use per job (default 16), and the timeout per job in seconds (default 14400 = 2 hr)
#	The arguments for the calls should be accessible from the directory where this script is launched
#	This script should be launched as (for example; change depending on what the jobs need): 
#		qsub -l ncpus=192,mem=768GB,jobfs=240GB,walltime=12:00:00,storage=gdata/nm31,wd \
#			-v calls_file="/rel/path/to/calls.txt",cores_per="___",timeout="____" \
#			launch_parallel.sh
######################


# Set PBS directives
#PBS -P nm31
#PBS -q normal
#PBS -v calls_file,cores_per,timeout


# set up a master log file
logfile="$(pwd)"/parallel_jobs_"${PBS_JOBID/\.gadi-pbs/}".log


# load the module for handling parallel tasks
module load nci-parallel/1.0.0a >> $logfile 2>&1


# set how many cores per run
# per node there are 192GB/48cpus = 4 GB per cpu
# for example: ncpus=192, mem=768GB, cores_per=16 (64 GB per job; 12 concurrent)
if [ -z $cores_per ]; then
	export cores_per=16
else
	export cores_per="$cores_per"
fi


# check args
if [ -z "$calls_file" ]; then
	echo -e "\nPlease specify a calls file\n"
	exit 1
else
	calls_file="$(readlink -f $calls_file)"
fi

if [ -z "$timeout" ]; then
	timeout=14400
else
	timeout="$timeout"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting multiple jobs at $(date)\n**********" >> $logfile 2>&1


# submit the file with all the script calls to nci-parallel
# set timeout to longer than each job is expected to take
# (e.g., for 12 hour walltime, with 12 concurrent, and ~153 samples, need each sample to take no longer than ~50 minutes = 3000)
mpirun -np $((PBS_NCPUS/cores_per)) --map-by ppr:$((48/cores_per)):node:PE=${cores_per} \
	nci-parallel --input-file "$calls_file" --timeout "$timeout" >> $logfile 2>&1


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished all parallel jobs at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
