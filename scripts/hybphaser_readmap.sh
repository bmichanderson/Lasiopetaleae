#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: May 2022
# Modified: July 2022, Nov 2022 (reworked, optimised, made more general), Jan 2023 (added more parameters; changed output)
# Description:	run read mapping and consensus generation for non-HybPiper results into folder structure expected by HybPhaser
#	It assumes that the outpute will be run through HybPhaser using the "intronerated" option
# Note: requires a Singularity container with all the dependencies and scripts from Lars
#	Results will be synced to a folder "HybPhaser/01_data" in the launch directory
#	Should be run via a launch script (independent calls to this script), but could also be run as:
#		qsub -v read1="_____",read2="______",multifasta="_____" hybphaser_readmap.sh
#	The multifasta should be named with the sample number prior to the first underscore,
#	and entries in it should take the form ">sample-locus"
######################


# Set PBS directives for if calling independently
#PBS -v read1,read2,multifasta


# set cores if needed
if [ -z $cores_per ]; then
	# this was launched independently from the launch script
	cores_per="$PBS_NCPUS"
fi
# set max mem based on 2GB per core
max_mem=$(echo "${cores_per}*2*85/100" | bc)g


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/hybphaser.sif"


# Set script to split output into individual locus files
locus_script="/g/data/nm31/ben_anderson/scripts/split_gene_ref.py"


# bbmap args
minid="0.95"		# minimum identity for read mapping
ambig="random"		# the way to treat ambiguous mapping; random = randomly assign read; best = first best location
pairlen="400"		# maximum allowed distance between paired reads; insert = pairlen + read1 + read2
pairedonly="t"		# whether to map paired reads only
trimdescrip="t"		# whether to trim the read descriptions (important for later processing)

# consensus args (from Lars)
min_cover="10"			# minimum depth for assessing whether there should be an ambiguity code
min_allele_freq="0.15"	# minimum allele frequency to report an ambiguity
min_allele_count="4"	# minimum count of alternative allele to report an ambiguity


# read and check args
if [ -n "$read1" ] && [ -n "$read2" ] && [ -n "$multifasta" ]; then
	read1="$(readlink -f $read1)"
	read2="$(readlink -f $read2)"
	multifasta="$(readlink -f $multifasta)"	
elif [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo -e "\nPlease specify input reads and multifasta as arguments 1, 2 and 3\n"
	exit 1
else
	read1="$(readlink -f $1)"
	read2="$(readlink -f $2)"
	multifasta="$(readlink -f $3)"
fi


# change to the job directory, then make a directory for working
# setting the sample name based on the multifasta name
cd $PBS_JOBFS
file_name=$(basename ${multifasta})
sample=${file_name/_*/}
mkdir temp_"$sample" && cd temp_"$sample"
mkdir -p HybPhaser/01_data/"$sample"/intronerated_consensus
mkdir -p HybPhaser/01_data/"$sample"/intronerated_contigs
mkdir -p HybPhaser/01_data/"$sample"/reads
logfile="$(readlink -f HybPhaser/01_data/${sample}/${sample}.log)"


# report job start (in the overall screen output)
echo -e "This job is going to map to $sample using $cores_per cores with max mem of $max_mem"


# load the singularity and python modules
module load singularity python3/3.10.0 >> $logfile 2>&1


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser read mapping at $(date)\n**********" >> $logfile 2>&1


# run read mapping (creating "map_sorted.bam")
cp "$multifasta" ref.fasta
singularity exec -H "$(pwd)" "$container" bbmap.sh ref=ref.fasta in1="$read1" in2="$read2" \
	outm=map.sam bs=bs.sh nodisk minid="$minid" ambiguous="$ambig" pairedonly="$pairedonly" \
	pairlen="$pairlen" trimreaddescriptions="$trimdescrip" threads="$cores_per" \
	-Xmx"$max_mem" >> $logfile 2>&1
singularity exec -H "$(pwd)" "$container" bash bs.sh >> $logfile 2>&1
rm bs.sh map.sam


# create a script for the container to run that generates the consensus sequence
# this is based on Lars' generate_consensus script
search_line="(DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= $min_allele_freq && (DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= $min_cover && (DP4[2]+DP4[3]) >= $min_allele_count"
awk_line='{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
echo "bcftools mpileup -I -Ov -f ref.fasta map_sorted.bam | bcftools call -mv -A -Oz -o loci.vcf.gz &&" >> run.sh
echo "bcftools index -f --threads $cores_per loci.vcf.gz &&" >> run.sh
echo "bcftools consensus -I -i \"$search_line\" -f ref.fasta loci.vcf.gz | awk '$awk_line' > output.fasta" >> run.sh
echo "echo \"\" >> output.fasta" >> run.sh


# run the script with the container
singularity exec -H "$(pwd)" "$container" bash run.sh >> $logfile 2>&1


# extract the mapped reads per sample and move to the output folder
# the files need to be named "{locus}_combined.fasta"
mkdir temp1 && cd temp1
## grab the names of the reference loci in the bam file
mv ../map_sorted.bam* .
singularity exec -H "$(pwd)" "$container" samtools idxstats map_sorted.bam | cut -f 1 | sed '/^\*/d' > loci_names.txt
## for each name, grab the reads
for reflocus in $(cat loci_names.txt)
do
	singularity exec -H "$(pwd)" "$container" samtools view -h -b map_sorted.bam "$reflocus" > temp.bam
	singularity exec -H "$(pwd)" "$container" reformat.sh in=temp.bam out=temp.fasta >> $logfile 2>&1
	singularity exec -H "$(pwd)" "$container" repair.sh in=temp.fasta out="${reflocus/*-/}".fasta >> $logfile 2>&1
	rm temp.fasta temp.bam
done
## move all to the output folder
for file in *.fasta
do
	mv "$file" ../HybPhaser/01_data/"$sample"/reads/"${file/.fasta/_combined.fasta}"
done
cd ..


# process the "ref.fasta" and "output.fasta" into single files per entry and move to the output fodlers
mkdir temp2 && cd temp2
python3 "$locus_script" ../ref.fasta >> $logfile 2>&1
for file in *.fasta
do
	mv "$file" ../HybPhaser/01_data/"$sample"/intronerated_contigs/"${file/.fasta/_intronerated.fasta}"
done
cd ..

mkdir temp3 && cd temp3
python3 "$locus_script" ../output.fasta >> $logfile 2>&1
for file in *.fasta
do
	mv "$file" ../HybPhaser/01_data/"$sample"/intronerated_consensus/"${file/.fasta/_intronerated.fasta}"
done
cd ..


# remove temporary files
rm ref.fasta output.fasta
rm -r temp1/ temp2/ temp3/


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser read mapping at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1


# move the results back to where the job was launched from  
rsync -rut HybPhaser $PBS_O_WORKDIR/
cd ..
rm -r temp_"$sample"/


# report job completion (in the overall screen output)
echo -e "\tThe mapping and consensus generation for $sample finished in "$duration_mins" minutes"
