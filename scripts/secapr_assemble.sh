#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: June 2022
# Modified: Oct 2022
# Description: run the SPAdes assembly step of the SECAPR pipeline (creating expected folder structure first)
# Note: the main input argument is a directory where each sample has a folder with two cleaned fastq.gz read files
#		The folders should be named by sample, and the read files should end in *1.fastq.gz and *2.fastq.gz
# Note: it requires a Singularity container with the SECAPR pipeline from Tobias Andermann
# Note: it will use a custom (SPAdes default) kmer range (21,33,55,77) rather than the SECAPR default (adds 99,127)
# Note: results will be created in the current directory in a folder called contigs; launch as:
#		qsub -v qc_dir="_____" secapr_assemble.sh
######################

# Set PBS directives

#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=192
#PBS -l mem=384GB
#PBS -l jobfs=400GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v qc_dir


# Set cores
ncores="$PBS_NCPUS"


# Set max total memory based on 90% of 2GB per core
tot_max_mem=$(echo "${ncores}*2*90/100" | bc)


# Determine how many samples can be assembled in parallel (if we use, e.g., 8 cores per)
cores_per="8"
instances=$(echo "${ncores}/${cores_per}" | bc)
mem_per=$(echo "${tot_max_mem}/${instances}" | bc)


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/secapr.sif"


# Check args
if [ -z "$qc_dir" ]; then
	echo -e "\nPlease specify a directory with folders of samples containing cleaned reads\n"
	exit 1
else
	qc_dir="$(readlink -f $qc_dir)"
fi


# Set the output directory (typically the current directory)
out_dir="$(pwd)"


# set up a log file
logfile="$(pwd)"/secapr_assemble_"${PBS_JOBID/\.gadi-pbs/}".log


# Load the singularity module
module load singularity >> $logfile 2>&1


# Change to the working area
cd $PBS_JOBFS


# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting SECAPR assembly at $(date)\n**********" >> $logfile 2>&1


# Create the folder structure expected by SECAPR
mkdir "cleaned_reads"
for sampledir in "$qc_dir"/*/; do
	sample="$(basename $sampledir)"
	mkdir cleaned_reads/"$sample"
	ln -s "$sampledir"*1.fastq.gz cleaned_reads/"$sample"/"$sample"_0_clean-READ1.fastq.gz
	ln -s "$sampledir"*2.fastq.gz cleaned_reads/"$sample"/"$sample"_0_clean-READ2.fastq.gz
done


# Run the assemblies
singularity exec -H "$(pwd)" "$container" secapr assemble_reads --input cleaned_reads --output contigs \
	--max_memory "$mem_per" --cores "$cores_per" --instances "$instances" --kmer 21,33,55,77 >> $logfile 2>&1


# Sync the results back to the original directory
rsync -rut contigs $PBS_O_WORKDIR/


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished SECAPR assembly at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
