#!/bin/bash

######################
# Author: Benjamin Anderson, with ideas from Theo Allnut and his script "hp2.py"
# Date: Apr 2023
# Description: run HybPiper jobs on Gadi
# Note: It should be launched as: 
#		qsub -v cores_per="__[defualt 8]",target_file="__",reads_dir="__",samples_file="___",
#			target_type="dna[default] or aa",intronerate="n or y",add_args="___" hybpiper.sh
######################


# Set PBS directives

#PBS -P nm31
#PBS -q normal
#PBS -l ncpus=192
#PBS -l mem=384GB
#PBS -l jobfs=400GB
#PBS -l walltime=24:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v cores_per,target_file,reads_dir,samples_file,target_type,intronerate,add_args


# set up a master log file
logfile="$(pwd)"/hybpiper_"${PBS_JOBID/\.gadi-pbs/}".log


# load the module for handling parallel tasks and the Python3 module
module load nci-parallel/1.0.0a python3/3.10.0 >> $logfile 2>&1


# export the user's home bin directory to PATH
export PATH="$PATH":"$HOME"/.local/bin


# set how many cores per task
# per node there are 192 GB / 48 cpus = 4 GB per cpu; but this is not necessarily needed (go with 2 GB per)
# if using ncores_per_task=8 (16 GB per job; 24 concurrent)
if [ -z $cores_per ]; then
	export ncores_per_task=8
else
	export ncores_per_task="$cores_per"
fi


# check args
if [ -z "$target_file" ]; then
	echo -e "\nPlease specify a target file\n" >> $logfile 2>&1
	exit 1
else
	target_file="$(readlink -f $target_file)"
fi
if [ -z "$reads_dir" ]; then
	echo -e "\nPlease specify a reads directory\n" >> $logfile 2>&1
	exit 1
else
	reads_dir="$(readlink -f $reads_dir)"
fi
if [ -z "$samples_file" ]; then
	echo -e "\nPlease specify a samples file\n" >> $logfile 2>&1
	exit 1
else
	samples_file="$(readlink -f $samples_file)"
fi
if [ -z "$target_type" ]; then
	target_type="dna"
elif [ "$target_type" = "aa" ] || [ "$target_type" = "dna" ]; then
	target_type="$target_type"
else
	echo -e "\nPlease specify target type as \"aa\" or \"dna\"\n" >> $logfile 2>&1
fi
if [ -z "$intronerate" ]; then
	intronerate="n"
elif [ "$intronerate" = "y" ] || [ "$intronerate" = "n" ]; then
	intronerate="$intronerate"
else
	echo -e "\nPlease specify intronerate as \"y\" or \"n\"\n" >> $logfile 2>&1
fi


# create the argument array and append additional arguments if present
arg_array=()
arg_array+=("--diamond")
arg_array+=("--merged")
if [ "$intronerate" = "y" ]; then
	arg_array+=("--run_intronerate")
fi
if [ -n "$add_args" ]; then
	IFS=" " read -r -a temp_array <<< "$add_args"
	for element in "${temp_array[@]}"
	do
		arg_array+=("$element")
	done
fi

# determine the form of the calls from the target type
if [ "$target_type" = "dna" ]; then
	targ_piece="-t_dna"
else
	targ_piece="-t_aa"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting HybPiper runs at $(date)\n**********" >> $logfile 2>&1


# make required directories
if [ ! -d logs ]; then
	mkdir -p logs
fi
if [ ! -d results ]; then
	mkdir -p results
fi


# create a bash script to run the individual HybPiper assemblies for each sample  
# the script will take arg1 = sample name
# the script will also run sequence collection, since we can't store files on Gadi
echo -e '#!/bin/bash' > assemble_script.sh
echo -e 'start="$(date +%s)"' >> assemble_script.sh
echo -e 'echo -e "This job is going to assemble sample $1 using ${ncores_per_task} cores"' >> assemble_script.sh
echo >> assemble_script.sh
echo -e 'mkdir results/"$1"' >> assemble_script.sh
echo -e 'cd "$PBS_JOBFS"' >> assemble_script.sh
echo >> assemble_script.sh
# assemble
echo -e 'hybpiper assemble' "$targ_piece" "$target_file" '-r' "$reads_dir"'/"$1"_* --prefix "$1"' \
	'--cpu "$ncores_per_task"' "${arg_array[@]}" '>>' "$PBS_O_WORKDIR"'/logs/"$1".out 2>>' \
	"$PBS_O_WORKDIR"'/logs/"$1".err' >> assemble_script.sh
echo >> assemble_script.sh
echo -e 'end="$(date +%s)"' >> assemble_script.sh
echo -e 'duration="$(( $end - $start ))"' >> assemble_script.sh
echo -e 'duration_mins=$(echo "scale=2; ${duration}/60" | bc)' >> assemble_script.sh
echo -e 'echo -e "\tThe assembly for sample $1 finished in "$duration_mins" minutes"' >> assemble_script.sh
echo -e 'echo "$1" >> "$1"_done.txt' >> assemble_script.sh
echo >> assemble_script.sh
# collect stats
echo -e 'hybpiper stats' "$targ_piece" "$target_file" 'gene "$1"_done.txt --seq_lengths_filename' \
	"$PBS_O_WORKDIR"'/results/"$1"/"$1"_lengths --stats_filename' "$PBS_O_WORKDIR"'/results/"$1"/"$1"_stats >>' \
	"$PBS_O_WORKDIR"'/logs/"$1".out 2>>' "$PBS_O_WORKDIR"'/logs/"$1".err' >> assemble_script.sh
echo >> assemble_script.sh
# collect sequences
echo -e 'hybpiper retrieve_sequences' "$targ_piece" "$target_file" 'dna --sample_names "$1"_done.txt --fasta_dir' \
	"$PBS_O_WORKDIR"'/results/"$1"/dna_seqs >>' "$PBS_O_WORKDIR"'/logs/"$1".out 2>>' \
	"$PBS_O_WORKDIR"'/logs/"$1".err' >> assemble_script.sh
if [ "$intronerate" = "y" ]; then
	echo >> assemble_script.sh
	echo -e 'hybpiper retrieve_sequences' "$targ_piece" "$target_file" 'supercontig --sample_names "$1"_done.txt --fasta_dir' \
		"$PBS_O_WORKDIR"'/results/"$1"/supercontig_seqs >>' "$PBS_O_WORKDIR"'/logs/"$1".out 2>>' \
		"$PBS_O_WORKDIR"'/logs/"$1".err' >> assemble_script.sh
fi
echo >> assemble_script.sh
# collect paralogs
echo -e 'hybpiper paralog_retriever' "$targ_piece" "$target_file" '"$1"_done.txt --fasta_dir_all' \
	"$PBS_O_WORKDIR"'/results/"$1"/paralogs_all --fasta_dir_no_chimeras' "$PBS_O_WORKDIR"'/results/"$1"/paralogs_no_chimeras' \
	'--paralog_report_filename' "$PBS_O_WORKDIR"'/results/"$1"/paralog_report >>' "$PBS_O_WORKDIR"'/logs/"$1".out 2>>' \
	"$PBS_O_WORKDIR"'/logs/"$1".err' >> assemble_script.sh
echo >> assemble_script.sh
# report completion
echo -e 'echo "$1" >>' "$PBS_O_WORKDIR"'/done.txt' >> assemble_script.sh
echo >> assemble_script.sh
# remove unnecessary files and compress results
echo -e 'mv "$1"/genes_with_seqs.txt "$1"/temp' >> assemble_script.sh
echo -e 'rm "$1"/*.*' >> assemble_script.sh
echo -e 'rm "$1"/*/*.fastq' >> assemble_script.sh
echo -e 'rm "$1"/*/"$1"/*.*' >> assemble_script.sh
echo -e 'rm -r "$1"/*/"$1"/sequences/FAA' >> assemble_script.sh
echo -e 'mv "$1"/temp "$1"/genes_with_seqs.txt' >> assemble_script.sh
echo -e 'tar -c --use-compress-program=pigz -f' "$PBS_O_WORKDIR"'/results/"$1"/"$1".tar.gz "$1"' >> assemble_script.sh
echo -e 'rm -r "$1"*' >> assemble_script.sh


# create a calls file
for sample in $(cat $samples_file); do
	echo -e "bash assemble_script.sh $sample" >> calls.txt
done


# launch the jobs
# set timeout to longer than each job is expected to take
# (e.g. for 48 hour walltime, with 24 concurrent, and 144 samples, need each sample to take no longer than 8 hours = 28,800)
mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((48/ncores_per_task)):node:PE=${ncores_per_task} \
	nci-parallel --input-file calls.txt --timeout 28800 >> $logfile 2>&1


# remove extraneous lines from the log files
sed -i '/ETA:/d' logs/*.err


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished all HybPiper runs at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
