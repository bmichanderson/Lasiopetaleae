#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: May–June 2023
# Description: run a GetOrganelle assembly of a plastome from an input SPAdes assembly
# Note: a Singularity container with the dependencies ("getorganelle.sif") is required
#	The script expects pre-downloaded database files in a folder (.GetOrganelle) in 
#	the launch directory; download like so (need to use copyq for access to the internet):
#		
#	container="/g/data/nm31/ben_anderson/singularity/getorganelle.sif"
#	singularity exec -H "$(pwd)" "$container" get_organelle_config.py --config-dir "$(pwd)"/.GetOrganelle \
#		--add embplant_pt,embplant_mt
#
# Note: the argument is a path to a folder with a SPAdes assembly (includes `scaffolds.paths` and `assembly_graph.fastg`)
#	The output will be directed by default to a folder in the launch directory called "getorganelle_out"; this
#	can be specified with the second argument `outdir`; launch as:
#		qsub -v folder="___",outdir="____" plastome_from_assembly.sh
#
# Note: if launching independently rather than in a batch, call as:
#	qsub -q normal \
#		-l ncpus=8,mem=32GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd \
#		-v folder="____",outdir="____" \
#		plastome_from_assembly.sh
######################


# Set PBS directives
#PBS -P nm31
#PBS -v folder,outdir


# Set cores based on whether this was launched independently or from a batch script
if [ -z $cores_per ]; then
	# this was launched independently from the launch script
	cores_per="$PBS_NCPUS"
fi
echo "This job is going to use $cores_per cores"


# Set container to use
container="/g/data/nm31/ben_anderson/singularity/getorganelle.sif"


# Set launch directory variable (where the database folder should be)
launchdir="$(pwd)"


# Check args
if [ -z "$folder" ]; then
	if [ -z "$1" ]; then
		echo -e "\nPlease specify a folder with a SPAdes assembly!\n"
		exit 1
	fi
	folder="$(readlink -f $1)"
else
	folder="$(readlink -f $folder)"
fi

if [ -z "$outdir" ]; then
	if [ -z "$2" ]; then
		outdir="$(pwd)"/getorganelle_out
	else
		outdir="$(readlink -f $2)"
	fi
else
	outdir="$(readlink -f $outdir)"
fi


# Check output directory
if [ ! -d "$outdir" ]; then
	mkdir -p "$outdir"
fi


# Set up a log file
logfile="$outdir"/plastome_assembly_"${PBS_JOBID/\.gadi-pbs/}".log


# Load the singularity module
module load singularity >> $logfile 2>&1


# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting plastome assembly at $(date)\n**********\n" >> $logfile 2>&1


# Change to the working area, create a temporary directory for this sample, and copy over the databases
cd $PBS_JOBFS
if [ ! -d $(basename $outdir) ]; then
	mkdir $(basename $outdir)
else
	echo -e "\nDirectory exists!\n"
fi

cd $(basename $outdir)

if [ -d "$launchdir"/.GetOrganelle ]; then
	rsync -rut "$launchdir"/.GetOrganelle .
	ls -al . >> $logfile 2>&1
	if [ -d "$(pwd)"/.GetOrganelle ]; then
		echo -e "\nThe folder .GetOrganelle should be present" >> $logfile 2>&1
	else
		echo -e "\nMissing the GetOrganelle databases! Something went wrong!" >> $logfile 2>&1
		exit 1
	fi
else
	echo -e "\nUnable to locate the GetOrganelle databases!" >> $logfile 2>&1
	echo -e "Please follow the instructions in the script for downloading before running the job" >> $logfile 2>&1
	exit 1
fi


# Run the assembly
# GetOrganelle
singularity exec -H "$(pwd)" "$container" get_organelle_from_assembly.py --config-dir "$(pwd)"/.GetOrganelle \
	-g "$folder"/assembly_graph.fastg -o getorganelle_out --spades-out-dir "$folder" \
	-F embplant_pt -t "$cores_per" > getorganelle.log 2>&1


# Sync the results to the output directory
rsync -rut getorganelle_out $outdir
rsync -rut *.log $outdir


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished plastome assembly at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
echo "This job finished"
