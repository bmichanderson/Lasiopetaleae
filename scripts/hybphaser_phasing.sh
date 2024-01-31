#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: Jan 2023
# Description: run HybPhaser phasing steps based on the clade association on Gadi
# Note: it requires a Singularity container with all the dependencies and scripts from Lars
# 		it requires a config.txt file in the out_dir specifying the required parameters (see documentation on github)  
# Note: it requires a file structure matching the config.txt file, and a csv list of accessions and references to phase to
#		Launch with args for the output directory (where the file structure is, e.g. has the folder 01_data present;
#		default is the launch dir)
#		qsub -v out_dir="_____" hybphaser_phasing.sh
# Note:	it will create folders "05_phasing" and "phased_reads" in the out_dir, so the config.txt file should match that
# Note:	set the number of threads and max memory in the config file to match this script (85% of memory in this script)
######################


# Set PBS directives

#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=24
#PBS -l mem=96GB
#PBS -l walltime=24:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v out_dir


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/hybphaser.sif"


# set up a master log file
logfile="$(pwd)"/hybphaser_phasing_"${PBS_JOBID/\.gadi-pbs/}".log


# Load the singularity module
module load singularity >> $logfile 2>&1


# Check args
if [ -z "$out_dir" ]; then
	out_dir="$(pwd)"
else
	out_dir="$(readlink -f $out_dir)"
fi

config_file="$out_dir"/config.txt
if [ ! -f $config_file ]; then
	echo -e "\nMissing configuration file (config.txt)\n"
	exit 1
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser phasing at $(date)\n**********" >> $logfile 2>&1


# Change to the outdirectory and run the R script to prepare the bbsplit script  
cd $out_dir
if [ ! -d "05_phasing" ]; then
	mkdir "05_phasing"
fi
singularity run -H "$(pwd)" "$container" 3a_prepare_phasing_script.R >> $logfile 2>&1


# Run the script that was created  
singularity exec -H "$(pwd)" "$container" bash "05_phasing"/run_bbsplit4phasing.sh


# Run the R script to process the output of the bbsplit runs  
singularity run -H "$(pwd)" "$container" 3b_collate_phasing_stats.R >> $logfile 2>&1


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser phasing at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
