#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: July 2022
# Modified: Oct 2022 (add intronerate option), Jan 2023 (added bbmap option)
# Description: run the first HybPhaser step on Gadi, using Lars' script in my singularity container
# Note: input will be a file with a list of samples, one per line, and the path to the HybPiper output (decompressed)
#		also, set whether this is from Theo's wrapper (needs decompression) [default] or not
#		also, set whether to run with the intronerate supercontigs or not [default]
#		also, set whether to use BBMap or not [default]
#		launch (where you want output folder) as: 
#		qsub -v samples="____",hp2_dir="____",decompress="y or n",intronerate="n or y",bbmap="n or y" hybphaser_generate.sh
#		the script will generate a folder called HybPhaser in the launch directory
######################

# Set PBS directives

#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=24
#PBS -l mem=192GB
#PBS -l jobfs=400GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v samples,hp2_dir,decompress,intronerate,bbmap


# define the container being used
container="/g/data/nm31/ben_anderson/singularity/hybphaser.sif"


# set the script for replacing NNNs when intronerate is set
nscript="/g/data/nm31/ben_anderson/scripts/fasta_Ns_sub.py"


# define the maximum number of samples to run at once (<= ncpus)
maxcontemp=24


# define the minimum identity for BBMap read mapping (if chosen)
minid=0.95


# set up a log file
logfile="$(pwd)"/hybphaser_generate_"${PBS_JOBID/\.gadi-pbs/}".log


# check args
if [ -z "$samples" ]; then
	echo -e "\nPlease specify a samples file\n"
	exit 1
else
	samples="$(readlink -f $samples)"
fi

if [ -z "$hp2_dir" ]; then
	echo -e "\nPlease specify a HybPiper results directory\n"
	exit 1
else
	hp2_dir="$(readlink -f $hp2_dir)"
fi

if [ -z "$decompress" ]; then
	decompress="y"
elif [ "$decompress" == "n" ] || [ "$decompress" == "y" ]; then
	decompress="$decompress"
else
	echo -e "\nPlease specify decompress as y or n\n"
	exit 1
fi

if [ -z "$intronerate" ]; then
	intronerate="n"
elif [ "$intronerate" == "n" ] || [ "$intronerate" == "y" ]; then
	intronerate="$intronerate"
else
	echo -e "\nPlease specify intronerate as n or y\n"
	exit 1
fi

if [ -z "$bbmap" ]; then
	bbmap="n"
elif [ "$bbmap" == "n" ] || [ "$bbmap" == "y" ]; then
	bbmap="$bbmap"
else
	echo -e "\nPlease specify bbmap as n or y\n"
	exit 1
fi

# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser consensus generation at $(date)\n**********" >> $logfile 2>&1


# load the singularity module
module load singularity python3/3.10.0 >> $logfile 2>&1


# Change to the working area
cd $PBS_JOBFS


# To make sure it doesn't go over size, split the samples into groups
split -l "$maxcontemp" "$samples" subset_
numbatches=$(find . -maxdepth 1 -name "subset_*" -printf '.' | wc -m)
echo -e "\nWill process $numbatches batches" >> $logfile 2>&1


# Move the HybPiper results to the working area
# If necessary, decompress the HybPiper results folders
index=1
if [ "$decompress" == "y" ]; then
	# Since there will be a lot to decompress, let's run in batches
	# Updated now to just run based on the split files (comment out old)
	#max_jobs="$PBS_NCPUS"
	for subset in subset_*
	do
		#cur_jobs=0
		echo -e "\n\nRunning decompression for batch $index" >> $logfile 2>&1
		for sample in $(cat $subset)
		do
			#((cur_jobs >= max_jobs)) && wait -n
			tar -xzf "$hp2_dir"/"$sample"/"$sample".tar.gz &
			#((++cur_jobs))
		done
		wait
		# Launch HybPhaser consensus generation
		echo -e "Running consensus generation for batch $index" >> $logfile 2>&1
		if [ "$intronerate" == "y" ]; then
			# first, we need to change the supercontigs to have 100 Ns instead of 10 when stitched
			echo -e "Substituting 100 Ns for 10 Ns in supercontig fasta files..." >> $logfile 2>&1
			python3 "$nscript" -m 100 */*/*/sequences/intron/*_supercontig.fasta
			# now, run the consensus generation
			echo -e "Generating consensus sequences...\n" >> $logfile 2>&1
			if [ "$bbmap" == "y" ]; then
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser -i -b -m "$minid" >> $logfile 2>&1
			else
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser -i >> $logfile 2>&1
			fi
		else
			echo >> $logfile 2>&1
			if [ "$bbmap" == "y" ]; then
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser -b -m "$minid" >> $logfile 2>&1
			else
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser >> $logfile 2>&1
			fi
		fi
		# manually remove mapping files and directories
		rm -r ./HybPhaser/*/*/mapping_files/
		# move the HybPhaser results to the working directory
		rsync -rut HybPhaser $PBS_O_WORKDIR/
		# remove the remaining files
		rm -r HybPhaser
		for sample in $(cat $subset)
		do
			rm -r "$sample"
		done
		((++index))
	done
else
	for subset in subset_*
	do
		for sample in $(cat $subset)
		do
			rsync -rt "$hp2_dir"/"$sample" .
		done
		# Launch HybPhaser consensus generation
		echo -e "\n\nRunning consensus generation for batch $index" >> $logfile 2>&1
		if [ "$intronerate" == "y" ]; then
			# first, we need to change the supercontigs to have 100 Ns instead of 10 when stitched
			echo -e "Substituting 100 Ns for 10 Ns in supercontig fasta files..." >> $logfile 2>&1
			python3 "$nscript" -m 100 */*/*/sequences/intron/*_supercontig.fasta
			# now, run the consensus generation
			echo -e "Generating consensus sequences...\n" >> $logfile 2>&1
			if [ "$bbmap" == "y" ]; then
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser -i -b -m "$minid" >> $logfile 2>&1
			else
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser -i >> $logfile 2>&1
			fi
		else
			echo >> $logfile 2>&1
			if [ "$bbmap" == "y" ]; then
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser -b -m "$minid" >> $logfile 2>&1
			else
				singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
					-t "$PBS_NCPUS" -o HybPhaser >> $logfile 2>&1
			fi
		fi
		# manually remove mapping files and directories
		rm -r ./HybPhaser/*/*/mapping_files/
		# move the HybPhaser results to the working directory
		rsync -rut HybPhaser $PBS_O_WORKDIR/
		# remove the remaining files
		rm -r HybPhaser
		for sample in $(cat $subset)
		do
			rm -r "$sample"
		done
		((++index))
	done
fi


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser consensus generation at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
