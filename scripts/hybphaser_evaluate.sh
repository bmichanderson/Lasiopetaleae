#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: May 2022
# Modified: July 2022, Oct 2022 (add intronerate option), Nov 2022 (add seqlist)
# Description: run HybPhaser evaluation of the results of a HybPiper run on Gadi
# Note: it requires a Singularity container with all the dependencies and scripts from Lars
# 		If wanting to customize the run, create a config.txt file as desired in the out_dir (then don't need the other args)
# Note: it also requires a file structure and consensus reads already present (from hybphaser_readmap.sh or hybphaser_generate.sh)
#		Launch with args for the output directory (where the file structure is, e.g. has the folder 01_data present),
#		the targets file (used for the HP2 run), and the namelist file (one sample ID per line)
#		also, specify whether the output to evaluate is intronerate or not [default]
#		also, specify whether to create sequence lists (if a custom config is provided for filtering)
#		qsub -v out_dir="_____",target_file="____",namelist="____",intronerate="n or y",seqlist="n or y" hybphaser_evaluate.sh
######################


# Set PBS directives

#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v out_dir,target_file,namelist,intronerate,seqlist


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/hybphaser.sif"


# set up a master log file
logfile="$(pwd)"/hybphaser_evaluate_"${PBS_JOBID/\.gadi-pbs/}".log


# load the singularity module
module load singularity >> $logfile 2>&1


# check args
if [ -z "$out_dir" ]; then
	out_dir="$(pwd)"
else
	out_dir="$(readlink -f $out_dir)"
fi

config_file="$out_dir"/config.txt
if [ ! -f $config_file ]; then
	if [ -z "$target_file" ]; then
		echo -e "\nPlease specify a target file\n"
		exit 1
	else
		target_file="$(readlink -f $target_file)"
	fi
	if [ -z "$namelist" ]; then
		echo -e "\nPlease specify a namelist file\n"
		exit 1
	else
		namelist="$(readlink -f $namelist)"
	fi
	if [ -z "$intronerate" ]; then
		intronerate="no"
	elif [ "$intronerate" == "n" ]; then
		intronerate="no"
	elif [ "$intronerate" == "y" ]; then
		intronerate="yes"
	else
		echo -e "\nPlease specify intronerate as n or y\n"
		exit 1
	fi
fi

if [ -z "$seqlist" ]; then
	seqlist="no"
elif [ "$seqlist" == "n" ]; then
	seqlist="no"
elif [ "$seqlist" == "y" ]; then
	seqlist="yes"
else
	echo -e "\nPlease specify seqlist as n or y\n"
	exit 1
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser evaluation at $(date)\n**********" >> $logfile 2>&1


# We need to create the config.txt file in the out_dir if it doesn't exist already
if [ ! -f $config_file ]; then
	echo "# General settings" > $config_file
	echo "path_to_output_folder = \"${out_dir}\"" >> $config_file
	echo "fasta_file_with_targets = \"${target_file}\"" >> $config_file
	echo "targets_file_format = \"DNA\"" >> $config_file
	echo "path_to_namelist = \"${namelist}\"" >> $config_file
	echo "intronerated_contig = \"$intronerate\"" >> $config_file
	echo "" >> $config_file

	echo "name_for_dataset_optimization_subset = \"\"" >> $config_file
	echo "" >> $config_file
	echo "# Missing data" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_loci_recovered = 0" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered = 0" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_samples_recovered = 0" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered = 0" >> $config_file
	echo "" >> $config_file

	echo "# Paralogs" >> $config_file
	echo "remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs = \"none\"" >> $config_file
	echo "file_with_putative_paralogs_to_remove_for_all_samples = \"\"" >> $config_file
	echo "remove_outlier_loci_for_each_sample = \"no\"" >> $config_file
	echo "" >> $config_file
else
	echo -e "\nDetected a config file, so will use it for the run\n" >> $logfile 2>&1
fi


# Change to the outdirectory and run the R script to count SNPs
cd $out_dir
singularity run -H "$(pwd)" "$container" 1a_count_snps.R >> $logfile 2>&1


# Run the R script to assess the dataset
singularity run -H "$(pwd)" "$container" 1b_assess_dataset.R >> $logfile 2>&1


# If requested, run the sequence list generation
if [ "$seqlist" == "yes" ]; then
	singularity run -H "$(pwd)" "$container" 1c_generate_sequence_lists.R >> $logfile 2>&1
fi


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser evaluation at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
