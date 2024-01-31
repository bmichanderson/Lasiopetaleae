#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: Jun 2023
# Description: run a SPAdes assembly using Illumina paired-end reads (or two sets of reads)
# Note: a Singularity container with SPAdes installed ("getorganelle.sif") is required
#	the argument is a string of semicolon-separated paths to up to two sets of fastq files (may be gzipped):
#		"/path/to/forward1.fq;/path/to/reverse1.fq[;/path/to/forward2.fq;/path/to/reverse2.fq]"
#	launch by itself as:
#		qsub -q normal \
#			-l ncpus=16,mem=64GB,jobfs=20GB,walltime=06:00:00,storage=gdata/nm31,wd \
#			-v reads="__;__" \
#			spades.sh
######################


# Set PBS directives
#PBS -P nm31
#PBS -v reads


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/getorganelle.sif"


# Set cores based on whether this was launched independently or from a batch script
if [ -z $cores_per ]; then
	# this was launched independently from the launch script
	cores_per="$PBS_NCPUS"
fi


# Set OMP variable for SPAdes multithreading
export OMP_NUM_THREADS="$cores_per"


# Set maximum memory for the job (ensure at least 4GB available per core)
max_mem=$(echo "${cores_per}*4*90/100" | bc)
echo "This job is going to use $cores_per cores with max mem of ${max_mem}GB"


# Check args
if [ -z "$reads" ]; then
	if [ -z "$1" ]; then
		echo -e "\nPlease specify paths to reads files!\n"
		exit 1
	fi
	reads="$1"
else
	reads="$reads"
fi

IFS=";" read -r -a temp_array <<< $reads
reads1forward=$(readlink -f "${temp_array[0]}")
reads1reverse=$(readlink -f "${temp_array[1]}")
if [ ${#temp_array[@]} -gt 2 ]; then	# if more than one read set
	oneset="False"
	reads2forward=$(readlink -f "${temp_array[2]}")
	reads2reverse=$(readlink -f "${temp_array[3]}")
else	# one set of read files
	oneset="True"
fi


# Change to the working directory, create a directory for output and start a log file
# it will be based on the first read file name before an underscore
cd $PBS_JOBFS
file_name=$(basename ${reads1forward})
prefix=${file_name/_*/}
if [ ! -d "$prefix" ]; then
	mkdir "$prefix"
fi
cd "$prefix"
logfile="$prefix".log


# Load the singularity module
module load singularity >> $logfile 2>&1


# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting SPAdes assembly at $(date)\n**********" >> $logfile 2>&1


# Run a SPAdes assembly differently depending on whether there are two sets of reads
if [ "$oneset" = "False" ]; then
	echo -e "\nDetected two sets of read files..." >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" spades.py --pe1-1 "$reads1forward" --pe1-2 "$reads1reverse" \
		--pe2-1 "$reads2forward" --pe2-2 "$reads2reverse" -o spades_out --only-assembler \
		-t "$cores_per" --memory "$max_mem" -k 21,45,65,85,105,125 > spades.log 2>&1
else
	echo -e "\nDetected one set of read files..." >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" spades.py -1 "$reads1forward" -2 "$reads1reverse" \
		-o spades_out --only-assembler -t "$cores_per" --memory "$max_mem" \
		-k 21,45,65,85,105,125 > spades.log 2>&1
fi


# Remove extra spades output, keeping the scaffolds and assembly graphs
mv "$(pwd)"/spades_out/scaffolds.fasta .
mv "$(pwd)"/spades_out/scaffolds.paths .
mv "$(pwd)"/spades_out/assembly_graph_with_scaffolds.gfa .
mv "$(pwd)"/spades_out/assembly_graph.fastg .
rm -r spades_out/


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished SPAdes assembly at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1


# Sync the results back to the launch directory
cd ..
rsync -rut "$prefix" $PBS_O_WORKDIR/

echo "This job finished"
