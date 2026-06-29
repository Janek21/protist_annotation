#!/bin/bash

#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

#SBATCH --job-name=ref_busco

#SBATCH --qos=normal
#SBATCH --time=90

#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

#SBATCH --array=0-57%5

#start
start_time=$(date +%s)
echo ">STARTING at $(date)"

allSpecies="$1" #list of species names(not taxid) as they appear in the folder

#How many lines to process per array task
LINES_PER_TASK=1

#Line range for this specific array ID
#If ID=0, START=1, END=4.
START_LINE=$(( SLURM_ARRAY_TASK_ID * LINES_PER_TASK + 1 ))
END_LINE=$(( (SLURM_ARRAY_TASK_ID + 1) * LINES_PER_TASK ))

echo "Processing lines $START_LINE through $END_LINE out of $(wc -l "$allSpecies")"

#Extract block of 4 SRR IDs and loop them
selected_specie=$(sed -n "${START_LINE},${END_LINE}p" $allSpecies)

echo "Species is $selected_specie"

#export so protein_busco can read
export SLURM_CPUS_PER_TASK

bash scripts/protein_busco.sh $selected_specie
echo "#################################################################################"
bash scripts/genome_busco.sh $selected_specie

rm agat_log*$selected_specie*

# Record memory usage (at the end of all 4 downloads)
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
# Check if the path exists to avoid errors on different cgroup versions
if [ -f "/sys/fs/cgroup$cgroup_dir/memory.peak" ]; then
	peak_mem=$(cat "/sys/fs/cgroup$cgroup_dir/memory.peak")
	peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}")
	echo ">Peak memory was $peak_mem_mb MegaBytes"
fi

#record end
elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60 )) minutes"
echo ">ENDING at $(date)"

