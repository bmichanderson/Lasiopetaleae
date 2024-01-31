#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: May 2022
# Modified: May 2023 (made deletion default), April 2023 (made Theo's script optional)
# Description: combine results from multiple HybPiper2 runs  
# Warning: the delete option will delete the original results after collating!
# NOTE: submit this from the directory where the folder "results" is:
#		qsub -v theo="n or y",intronerate="n or y",delete="n or y" combine_results.sh
######################


# Set PBS directives
#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=8GB
#PBS -l walltime=08:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v theo,intronerate,delete


# set locations of my scripts
len_script="/g/data/nm31/ben_anderson/scripts/hybpiper_loci_lengths.py"
comb_script="/g/data/nm31/ben_anderson/scripts/fasta_combining.py"


# Check args
if [ -z "$theo" ]; then
	theo="n"
elif [ "$theo" = "y" ] || [ "$theo" = "n" ]; then
	theo="$theo"
else
	echo -e "Please specify argument \"theo\" as \"y\" or \"n\""
	exit 1
fi

if [ -z "$intronerate" ]; then
	intronerate="n"
elif [ "$intronerate" = "y" ] || [ "$intronerate" = "n" ]; then
	intronerate="$intronerate"
else
	echo -e "Please specify argument \"intronerate\" as \"y\" or \"n\""
	exit 1
fi

if [ -z "$delete" ]; then
	delete="n"
elif [ "$delete" = "y" ] || [ "$delete" = "n" ]; then
	delete="$delete"
else
	echo -e "Please specify argument \"delete\" as \"y\" or \"n\""
	exit 1
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nCombining results at $(date)\n**********"


# load the python module
module load python3/3.10.0


# run Theo's script if requested; else run my implementation
if [ "$theo" = "y" ]; then
	join_hp_results_para.py results/ $PBS_NCPUS
else
	# change into the results directory
	cd results

	# combine stats files into one
	ls -d */ | sed 's/\///' > temp_samples.txt
	sample=$(head -n 1 temp_samples.txt)
	head -n 1 "$sample"/*lengths.tsv | cut -f 2- | tr '\t' '\n' | sort -n > temp_loci.txt
	head -n 1 "$sample"/*stats.tsv > combined_stats.tsv
	for samp in $(cat temp_samples.txt); do 
		tail -1 "$samp"/*stats.tsv >> combined_stats.tsv
	done

	# combine lengths files into one
	# NOTE: the gene lengths are in a different order for each sample
	# NOTE: the mean lengths are incorrect (different for each)
	python3 "$len_script" */*lengths.tsv

	# combine the fasta files into one per locus
	if [ ! -d dna_seqs ]; then
		mkdir dna_seqs
	fi
	if [ ! -d paralog_seqs ]; then
		mkdir paralog_seqs
	fi
	if [ "$intronerate" = "y" ]; then
		if [ ! -d supercontig_seqs ]; then
			mkdir supercontig_seqs
		fi
	fi

	# create a function for combining files
	combine_fastas () {
		loc="$1"
		python3 "$comb_script" -o "$loc"_temp.fasta */dna_seqs/"$loc".FNA
		mv "$loc"_temp.fasta dna_seqs/"$loc".fasta

		python3 "$comb_script" -o "$loc"_temp2.fasta */paralogs_all/"$loc"_paralogs_all.fasta
		mv "$loc"_temp2.fasta paralog_seqs/"$loc".fasta

		if [ "$intronerate" = "y" ]; then
			python3 "$comb_script" -o "$loc"_temp3.fasta */supercontig_seqs/"$loc"_supercontig.fasta
			mv "$loc"_temp3.fasta supercontig_seqs/"$loc".fasta
		fi
	}

	# to parallelise, split the loci into groups
	split -l "$PBS_NCPUS" temp_loci.txt temp_subset_
	for subfile in temp_subset_*; do
		for loc in $(cat $subfile); do
			combine_fastas "$loc" &
		done
		wait
	done

	# remove the individual sequences if requested
	if [ "$delete" = "y" ]; then
		rm -r */dna_seqs/
		rm -r */paralogs_all/
		rm -r */paralogs_no_chimeras/
		if [ "$intronerate" = "y" ]; then
			rm -r */supercontig_seqs/
		fi
		rm temp*
	fi
fi


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished combining results after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours"
