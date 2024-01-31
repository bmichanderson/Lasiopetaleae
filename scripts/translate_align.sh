#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: April 2023
# Description: translate multifastas that are in-frame, then align with MAFFT and convert back to nucleotide  
# Note:	it requires a Singularity container (phylo.sif) with the dependencies (MAFFT and pal2nal),
#		an argument specifying a text file of samples to drop from each fasta (one sample per line; optional),
#		and a folder with the input multifastas
#		It assumes the input multifastas are named *_name.fasta (so, will keep the last field separated by "_"
#		or everything if no underscore)
# Note:	submit from a directory where you want the resulting alignments as:
#		qsub -v dropfile="_____",fasta_dir="_____" translate_align.sh
######################


# Set PBS directives
#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v dropfile,fasta_dir


# set the container and scripts
container="/g/data/nm31/ben_anderson/singularity/phylo.sif"
remove_script="/g/data/nm31/ben_anderson/scripts/remove_fastas.py"
translate_script="/g/data/nm31/ben_anderson/scripts/translate.py"


# set up a master log file
logfile="$(pwd)"/translate_align_"${PBS_JOBID/gadi-pbs/}"log


# load the singularity and python3 modules
module load singularity python3/3.10.0 >> $logfile 2>&1


# check args
if [ -z "$fasta_dir" ]; then
	echo -e "\nPlease specify where the fasta files are\n"
	exit 1
else
	fasta_dir="$(readlink -f $fasta_dir)"
fi

if [ -z "$dropfile" ]; then
	drop="n"
else
	dropfile="$(readlink -f $dropfile)"
	drop="y"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting translation and alignment at $(date)\n**********" >> $logfile 2>&1


# if samples need to be dropped, do so, otherwise copy files to current directory
if [ "$drop" = "y" ]; then
	for fasta_file in "$fasta_dir"/*.fasta; do
		python3 "$remove_script" -f "$dropfile" "$fasta_file" >> $logfile 2>&1
		for file in mod*.fasta; do 
			newname=$(echo $file | awk -F_ '{print $NF}'); mv $file $newname
		done
	done
else
	for fasta_file in "$fasta_dir"/*.fasta; do
		newname=$(basename $fasta_file | awk -F_ '{print $NF}')
		cp "$fasta_file" ./"$newname"
	done
fi


# translate the files and align them
python3 "$translate_script" *.fasta >> $logfile 2>&1
sed -i 's/\*/X/g' *_prot.fasta		# convert stop codons to "X"
for translation in *_prot.fasta; do
	singularity exec "$container" mafft --auto --thread "$PBS_NCPUS" "$translation" > "${translation/\.fasta/\.aln}"
done


# convert back to nucleotide
for alignment in *_prot.aln; do
	singularity exec "$container" pal2nal "$alignment" "${alignment/_prot\.aln/\.fasta}" \
		-output fasta > "${alignment/prot\.aln/exon_aligned\.fasta}"
done


# rename final files to locus names and remove intermediate files
for file in *exon_aligned.fasta; do
	mv "$file" "${file/_exon_aligned/}"
	rm "${file/exon_aligned\.fasta/prot}"*
done


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished translation and alignment run at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
