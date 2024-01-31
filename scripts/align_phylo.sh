#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: May 2022
# Modified: Aug 2022, Nov 2022 (gappyout), Jan 2023 (more control over which analyses to run), Mar 2023 (custom cleaning script),
#		Mar 2023 (add locus bootstrapping and likelihood mapping), Jun 2023 (corrections and clean up)
# Description: run alignment clean up and phylogenetic analysis for a directory of fasta gene alignments
# Note:	it requires a Singularity container with all the dependencies/programs,
#		an argument indicating whether or not to re-align the input [default re-align], an argument for
#		whether to clean the alignments before analysis [default clean], an argument for process control
#		indicating which analyses to run: a=all c=concat l=loci [default "a"; other options: "c" "l"], an argument
#		indicating whether to run likelihood mapping [default do not], and an argument for whether to remove low support
#		branches in locus trees and turn them into polytomies [default do not]
# Note:	submit from a directory where you want the results as:
#		qsub -v align_dir="_____",realign="y or n",clean="y or n",analysis="a, c or l",likemap="n or y",poly="n or y" align_phylo.sh
######################


# Set PBS directives
#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=48GB
#PBS -l walltime=24:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v align_dir,realign,clean,analysis,likemap,poly


# set the max memory for java for Astral
max_mem="24G"


# set the percentage threshold for removing low support branches in locus trees for concordance factor and ASTRAL analyses
polytomy="50"


# set the container and cleaning script
container="/g/data/nm31/ben_anderson/singularity/phylo.sif"
clean_script="/g/data/nm31/ben_anderson/scripts/clean_alignment.py"


# set up a master log file
logfile=align_phylo_"${PBS_JOBID/gadi-pbs/}"log
logfile="$(readlink -f $logfile)"


# load the singularity and python3 modules
module load singularity python3/3.10.0 >> $logfile 2>&1


# check args
if [ -z "$align_dir" ]; then
	echo -e "\nPlease specify where the alignments are\n"
	exit 1
else
	align_dir="$(readlink -f $align_dir)"
fi

if [ -z "$realign" ]; then
	realign="y"
elif [ "$realign" == "n" ] || [ "$realign" == "y" ]; then
	realign="$realign"
else
	echo -e "\nPlease specify realign as y or n\n"
	exit 1
fi

if [ -z "$clean" ]; then
	clean="y"
elif [ "$clean" == "n" ] || [ "$clean" == "y" ]; then
	clean="$clean"
else
	echo -e "\nPlease specify clean as y or n\n"
	exit 1
fi

if [ -z "$analysis" ]; then
	analysis="a"
elif [ "$analysis" == "a" ] || [ "$analysis" == "c" ] || [ "$analysis" == "l" ]; then
	analysis="$analysis"
else
	echo -e "\nPlease specify analysis as a, c or l\n"
	exit 1
fi

if [ -z "$likemap" ]; then
	likemap="n"
elif [ "$likemap" == "n" ] || [ "$likemap" == "y" ]; then
	likemap="$likemap"
else
	echo -e "\nPlease specify likemap as y or n\n"
	exit 1
fi

if [ -z "$poly" ]; then
	poly="n"
elif [ "$poly" == "n" ] || [ "$poly" == "y" ]; then
	poly="$poly"
else
	echo -e "\nPlease specify poly as y or n\n"
	exit 1
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting alignment clean and phylogenetics run at $(date)\n**********" >> $logfile 2>&1


# make a directory for the alignments to be put in and to be input to phylogenetic analysis
in_dir="in_align"
if [ ! -d "$in_dir" ]; then
	mkdir -p "$in_dir"
fi


# re-align (if not specified not to)
if [ "$realign" == "y" ]; then
	echo -e "\n****\tRunning initial re-alignments" >> $logfile 2>&1
	align_start="$(date +%s)"
	# make a temp directory and re-align
	if [ ! -d temp_align ]; then
		mkdir temp_align
	fi

	for align in "$align_dir"/*.fasta
	do
		align_name="$(basename $align)"
		singularity exec -H "$(pwd)" "$container" mafft --quiet --auto --thread "$PBS_NCPUS" "$align" > temp_align/"$align_name"
	done
	align_end="$(date +%s)"
	duration="$(( $align_end - $align_start ))"
	echo -e "Finished re-alignments in $duration seconds\n" >> $logfile 2>&1
	align_dir=temp_align
else
	echo -e "\n****\tSkipping re-alignment\n" >> $logfile 2>&1
fi


# Clean the alignments (if not specified not to)  
# This used to be using Trimal with the "gappyout" automated mode (commented out; change if desired)  
# Now using my own script to remove positions with >50% gap/N and sequences with >75% gap/N
if [ "$clean" == "y" ]; then
	echo -e "\n****\tCleaning alignments\n" >> $logfile 2>&1
	clean_start="$(date +%s)"
	for align in "$align_dir"/*.fasta
	do
		align_name="$(basename $align)"
		echo -e "Cleaning alignment $align_name" >> $logfile 2>&1
#		singularity exec -H "$(pwd)" "$container" trimal -in "$align" -out "$in_dir"/"${align_name/.fasta/_clean.fasta}" \
#			-gappyout >> $logfile 2>&1
		ln -s "$align" .
		python3 "$clean_script" -p 50 -s 75 "$align_name" >> $logfile 2>&1
		mv "${align_name/.fasta/_clean.fasta}" "$in_dir"
		rm "$align_name"
	done
	clean_end="$(date +%s)"
	duration="$(( $clean_end - $clean_start ))"
	echo -e "\n****\tFinished cleaning alignments in $duration seconds\n" >> $logfile 2>&1
else
	echo -e "\n****\tSkipping cleaning step\n" >> $logfile 2>&1
	rsync -rut "$align_dir"/*.fasta "$in_dir"
fi


# remove the temporary alignments, if present
if [ -d temp_align ]; then
	rm -r temp_align
fi


# Run a concatenation tree using IQTREE (in the background if also running loci)
if [ "$analysis" == "a" ]; then
	# set the amount of threads to give to concatenation, and the remainder to loci
	threads_concat=8
	threads_loci=$((PBS_NCPUS - threads_concat))

	# run the analysis in the background so that loci can be started
	echo -e "\n****\tRunning concatenation tree in the background (using maximum $threads_concat CPUs)\n" >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" iqtree --threads-max "$threads_concat" -T AUTO \
		--ufboot 1000 --sampling GENESITE -m MFP --prefix concat -p "$in_dir" > /dev/null 2>&1  &
elif [ "$analysis" == "c" ]; then
	# run the analysis with all resources (not running loci)
	echo -e "\n****\tRunning concatenation tree (using maximum $PBS_NCPUS CPUs)\n" >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" iqtree --threads-max "$PBS_NCPUS" -T AUTO \
		--ufboot 1000 --sampling GENESITE -m MFP --prefix concat -p "$in_dir" > /dev/null 2>&1
else
	echo -e "\n****\tNot running concatenation tree\n" >> $logfile 2>&1
fi


# Run trees for each alignment
if [ "$analysis" == "a" ]; then
	# run the loci analyses while concat is running in the background
	echo -e "\n****\tStarting ML runs of all loci alignments (using maximum $threads_loci CPUs)\n" >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" iqtree --threads-max "$threads_loci" -T AUTO \
		--ufboot 1000 -m MFP --prefix loci -S "$in_dir" > /dev/null 2>&1
elif [ "$analysis" == "l" ]; then
	# run the loci analyses with all resources
	echo -e "\n****\tStarting ML runs of all loci alignments (using maximum $PBS_NCPUS CPUs)\n" >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" iqtree --threads-max "$PBS_NCPUS" -T AUTO \
		--ufboot 1000 -m MFP --prefix loci -S "$in_dir" > /dev/null 2>&1
else
	echo -e "\n****\tSkipping loci runs\n" >> $logfile 2>&1
fi


# Run likelihood mapping if requested
if [ "$likemap" == "y" ]; then
	mkdir likemap
	rsync -rut "$in_dir"/*.fasta likemap/
	cd likemap
	if [ "$analysis" == "a" ]; then
		# run the mapping while concat is still running in the background
		echo -e "\n****\tStarting likelihood mapping for loci alignments (using maximum $threads_loci CPUs)\n" >> $logfile 2>&1
		max_jobs=$((threads_loci/2))
	else
		# run the mapping with all resources
		echo -e "\n****\tStarting likelihood mapping for loci alignments (using maximum $PBS_NCPUS CPUs)\n" >> $logfile 2>&1
		max_jobs=$((PBS_NCPUS/2))
	fi
	curr_jobs=0
	for alignment in *.fasta
	do
		((curr_jobs >= max_jobs)) && wait -n
		singularity exec -H "$(pwd)" "$container" iqtree -T 2 -s "$alignment" \
			--lmap 10000 --quartetlh -n 0 --prefix lm_"${alignment/.fasta/}" > /dev/null 2>&1 &
		((++curr_jobs))
	done
	wait
	rm *.fasta *.phy *.log *.gz *.treefile *.eps *.svg
	cd ..
else
	echo -e "\n****\tSkipping likelihood mapping\n" >> $logfile 2>&1
fi


# wait for concatenation to finish if running in the background
if [ "$analysis" == "a" ]; then
	wait
	echo -e "...finished background job for the concatenation tree" >> $logfile 2>&1
fi


# Reduce low-support branches to polytomies, if requested
if [ "$poly" == "y" ]; then
	if [ -f "loci.treefile" ]; then
		echo -e "\n****\tReducing low supported branches (<$polytomy%) to polytomies" >> $logfile 2>&1
		singularity exec -H "$(pwd)" "$container" nw_ed loci.treefile "i & b < $polytomy" o > loci_poly.tre
		mv loci_poly.tre loci.treefile
	else
		echo -e "\n****\tNo loci tree file found; cannot reduce low-supported branches to polytomies" >> $logfile 2>&1
	fi
fi


# Calculate concordance factors on the concatenation tree if the relevant files exist
if [ -f "concat.treefile" ]; then
	if [ -f "loci.treefile" ]; then
		echo -e "\n****\tCalculating gene and site concordance factors" >> $logfile 2>&1
		singularity exec -H "$(pwd)" "$container" iqtree -T "$PBS_NCPUS" \
			-t concat.treefile --gcf loci.treefile --prefix concord_gcf \
			> /dev/null 2>&1
		singularity exec -H "$(pwd)" "$container" iqtree -T "$PBS_NCPUS" \
			-t concat.treefile --scfl 10000 --prefix concord_scf \
			-p "$in_dir" > /dev/null 2>&1
	else
		echo -e "\n****\tCalculating site concordance factors" >> $logfile 2>&1
		singularity exec -H "$(pwd)" "$container" iqtree -T "$PBS_NCPUS" \
			-t concat.treefile --scfl 10000 --prefix concord_scf \
			-p "$in_dir" > /dev/null 2>&1
	fi
fi


# Run an ASTRAL tree using the gene trees (if they exist)
if [ -f "loci.treefile" ]; then
	echo -e "\n****\tRunning ASTRAL analyses" >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" java -Xmx"$max_mem" -jar /Astral/astral.5.7.1.jar \
		-i loci.treefile -o astral.tre --branch-annotate 2 > astral_screenlog.txt 2>&1
	# add a polytomy test that will put p-values on branches
	singularity exec -H "$(pwd)" "$container" java -Xmx"$max_mem" -jar /Astral/astral.5.7.1.jar \
		-i loci.treefile -o astral_poly.tre --branch-annotate 10 >> astral_screenlog.txt 2>&1
else
	echo -e "\n****\tNo files found for ASTRAL analysis" >> $logfile 2>&1
fi


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished alignment clean and phylogenetics run at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
