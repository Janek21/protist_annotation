#!/bin/bash

#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#SBATCH --job-name=ref_busco

#SBATCH --qos=normal
#SBATCH --time=90

#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

#start
start_time=$(date +%s)
echo ">STARTING at $(date)"

#selected specie
selected_specie="$1"

echo "Species is $selected_specie"

#export so protein_busco can read
export SLURM_CPUS_PER_TASK

bash protein_busco.sh $selected_specie
echo "#################################################################################"
bash genome_busco.sh $selected_specie

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

