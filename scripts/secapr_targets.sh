#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: June 2022
# Description: run the search for target contigs step of the SECAPR pipeline
# Note: the arguments are a directory where each sample has a fasta file with contigs, and a targets fasta file
# Note: it requires a Singularity container with the SECAPR pipeline from Tobias Andermann
# Note: results will be created in the current directory in folders target_contigs and plots; launch as:
#		qsub -v contigs_dir="_____",targets_file="____" secapr_targets.sh
######################

# Set PBS directives

#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=8GB
#PBS -l walltime=02:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v contigs_dir,targets_file


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/secapr.sif"


# Check args
if [ -z "$contigs_dir" ]; then
	echo -e "\nPlease specify a directory with contig files\n"
	exit 1
else
	contigs_dir="$(readlink -f $contigs_dir)"
fi
if [ -z "$targets_file" ]; then
	echo -e "\nPlease specify a file with target fasta sequences, one per target\n"
	exit 1
else
	targets_file="$(readlink -f $targets_file)"
fi


# Make a log file
logfile="secapr_targets.log"


# Load the singularity module
module load singularity >> $logfile 2>&1


# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting SECAPR target search at $(date)\n**********" >> $logfile 2>&1


# Run the target search
singularity exec -H "$(pwd)" "$container" secapr find_target_contigs --contigs "$contigs_dir" \
	--reference "$targets_file" --target_length 100 --min_identity 85 \
	--output target_contigs --keep_paralogs >> $logfile 2>&1


# Plot a summary of recovered targets per sample
singularity exec -H "$(pwd)" "$container" secapr plot_sequence_yield --extracted_contigs target_contigs \
	--output plots >> $logfile 2>&1


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished SECAPR target search at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
