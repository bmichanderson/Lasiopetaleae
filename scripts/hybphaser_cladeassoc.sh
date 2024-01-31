#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: Jan 2023
# Modified: Jun 2023
# Description: run HybPhaser clade association steps on output from the evaluation steps on Gadi
# Note: it requires a Singularity container with all the dependencies and scripts from Lars
#	it requires a config.txt file in the out_dir specifying the required parameters (see documentation on github)
#	it requires the file structure produced from earlier steps of HybPhaser
#	Launch with args for the output directory (e.g. has the folder 01_data present; default = launch directory),
#	and an optional samples file listing the sample names (e.g. done.txt from HybPiper)
#	qsub -v out_dir="_____",samples="____" hybphaser_cladeassoc.sh
# Note:	it will create a folder "04_clade_association" in the out_dir, so the config.txt file should match that
# Note:	set the number of threads and max memory in the config.txt file to match this script (85% of memory in this script)
######################


# Set PBS directives
#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=24
#PBS -l mem=96GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v out_dir,samples


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/hybphaser.sif"


# Set up a master log file
logfile="$(pwd)"/hybphaser_cladeassoc_"${PBS_JOBID/\.gadi-pbs/}".log


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

if [ -z "$samples" ]; then
	sample_arg="True"
else
	sample_arg="False"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser clade association at $(date)\n**********" >> $logfile 2>&1


# Change to the outdirectory and run the script to gather mapped reads
cd $out_dir
if [ "$sample_arg" = "True" ]; then
	singularity exec -H "$(pwd)" "$container" 2_extract_mapped_reads.sh -n "$(readlink -f $sample_arg)" >> $logfile 2>&1	
else
	singularity exec -H "$(pwd)" "$container" 2_extract_mapped_reads.sh >> $logfile 2>&1
fi


# Run the R script to prepare the bbsplit script
if [ ! -d "04_clade_association" ]; then
	mkdir "04_clade_association"
fi
singularity run -H "$(pwd)" "$container" 2a_prepare_bbsplit_script.R >> $logfile 2>&1


# Run the script that was created
singularity exec -H "$(pwd)" "$container" bash "04_clade_association"/run_bbsplit4clade_association.sh


# Run the R script to process the output of the bbsplit runs
singularity run -H "$(pwd)" "$container" 2b_collate_bbsplit_results.R >> $logfile 2>&1


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser clade association at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
