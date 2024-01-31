#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: July 2023
# Description: run mapping followed by a SPAdes assembly of the mapped Illumina paired-end reads (or two sets of reads)
# Note: a Singularity container with SPAdes installed ("getorganelle.sif") is required
#	the arguments are a string of semicolon-separated paths to up to two sets of fastq files (may be gzipped):
#		"/path/to/forward1.fq;/path/to/reverse1.fq[;/path/to/forward2.fq;/path/to/reverse2.fq]"
#	a reference fasta to map to, and the minimum identity for mapping (default: 0.76)
#	launch by itself as, e.g.:
#		qsub -q normal \
#			-l ncpus=8,mem=16GB,jobfs=50GB,walltime=12:00:00,storage=gdata/nm31,wd \
#			-v reads="__;__",ref="___",ident="___" \
#			map_assemble.sh
######################


# Set PBS directives
#PBS -P nm31
#PBS -v reads,ref,ident


# define where bbmap is
bbmap_loc="/g/data/nm31/bin/bbmap"


# bbmap args
ambig="all"		# the way to treat ambiguous mapping; random = randomly assign read; best = first best location; all = retain all top-scoring sites


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/getorganelle.sif"


# Set cores based on whether this was launched independently or from a batch script
if [ -z $cores_per ]; then
	# this was launched independently from the launch script
	cores_per="$PBS_NCPUS"
fi


# Set max memory to specify (limit to 85% when passing to bbmap scripts; ensure at least 2GB available per core)
max_mem=$(echo "${cores_per}*2*85/100" | bc)


# Set OMP variable for SPAdes multithreading
export OMP_NUM_THREADS="$cores_per"


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

if [ -z "$ref" ]; then
	if [ -z "$2" ]; then
		echo -e "\nPlease specify reference fasta file!\n"
		exit 1
	fi
	ref="$(readlink -f $2)"
else
	ref="$(readlink -f $ref)"
fi

if [ -z "$ident" ]; then
	if [ -z "$3" ]; then
		echo -e "\nMinimum identity for read mapping not specified; using 0.76\n"
		minid=0.76
	fi
	minid="$3"
else
	minid="$ident"
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


# Report job starting
echo "This job (sample $prefix) is going to use $cores_per cores with max mem of ${max_mem}GB"


# Load the java and singularity modules
module load java/jdk-8.40 singularity >> $logfile 2>&1


# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting mapping and SPAdes assembly at $(date)\n**********" >> $logfile 2>&1


# Run mapping and assembly differently depending on if there are one or two sets of reads
if [ "$oneset" = "False" ]; then
	echo -e "\nDetected two sets of read files..." >> $logfile 2>&1
	"$bbmap_loc"/bbmap.sh ref="$ref" in1="$reads1forward" in2="$reads1reverse" outm=assemble1_R#.fastq.gz \
	nodisk minid="$minid" ambiguous="$ambig" threads="$cores_per" -Xmx"$max_mem"g >> $logfile 2>&1
	"$bbmap_loc"/bbmap.sh ref="$ref" in1="$reads1forward" in2="$reads1reverse" outm=map1.sam bs=bs1.sh \
	nodisk minid="$minid" ambiguous="$ambig" threads="$cores_per" -Xmx"$max_mem"g >> $logfile 2>&1
	./bs1.sh >> $logfile 2>&1

	"$bbmap_loc"/bbmap.sh ref="$ref" in1="$reads2forward" in2="$reads2reverse" outm=assemble2_R#.fastq.gz \
	nodisk minid="$minid" ambiguous="$ambig" threads="$cores_per" -Xmx"$max_mem"g >> $logfile 2>&1
	"$bbmap_loc"/bbmap.sh ref="$ref" in1="$reads2forward" in2="$reads2reverse" outm=map2.sam bs=bs2.sh \
	nodisk minid="$minid" ambiguous="$ambig" threads="$cores_per" -Xmx"$max_mem"g >> $logfile 2>&1
	./bs2.sh >> $logfile 2>&1

	rm bs*.sh map*.sam

	# assemble
	singularity exec -H "$(pwd)" "$container" spades.py --pe1-1 assemble1_R1.fastq.gz --pe1-2 assemble1_R2.fastq.gz \
		--pe2-1 assemble2_R1.fastq.gz --pe2-2 assemble2_R2.fastq.gz -o spades_out --only-assembler \
		-t "$cores_per" --memory "$max_mem" -k 21,45,65,85,105,125 > spades.log 2>&1
else
	echo -e "\nDetected one set of read files..." >> $logfile 2>&1
	"$bbmap_loc"/bbmap.sh ref="$ref" in1="$reads1forward" in2="$reads1reverse" outm=assemble1_R#.fastq.gz \
	nodisk minid="$minid" ambiguous="$ambig" threads="$cores_per" -Xmx"$max_mem"g >> $logfile 2>&1
	"$bbmap_loc"/bbmap.sh ref="$ref" in1="$reads1forward" in2="$reads1reverse" outm=map1.sam bs=bs1.sh \
	nodisk minid="$minid" ambiguous="$ambig" threads="$cores_per" -Xmx"$max_mem"g >> $logfile 2>&1
	./bs1.sh >> $logfile 2>&1

	rm bs*.sh map*.sam

	# assemble
	singularity exec -H "$(pwd)" "$container" spades.py -1 assemble1_R1.fastq.gz -2 assemble1_R2.fastq.gz \
		-o spades_out --only-assembler -t "$cores_per" --memory "$max_mem" \
		-k 21,45,65,85,105,125 > spades.log 2>&1
fi


# Remove extra spades output, keeping the scaffolds and assembly graphs
mv "$(pwd)"/spades_out/scaffolds.fasta .
mv "$(pwd)"/spades_out/scaffolds.paths .
mv "$(pwd)"/spades_out/assembly_graph_with_scaffolds.gfa .
mv "$(pwd)"/spades_out/assembly_graph.fastg .
rm -r spades_out/


# Remove the mapped reads
rm assemble*.fastq.gz


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished mapping and SPAdes assembly at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1


# Sync the results back to the launch directory
cd ..
rsync -rut "$prefix" $PBS_O_WORKDIR/


# Report job finishing
echo "This job (sample $prefix) finished"
