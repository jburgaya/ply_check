#!/bin/bash
#SBATCH --job-name plyCheck
#SBATCH --output=/storage/homefs/jb25u082/logs/plyCheck_Project1000.out
#SBATCH --error=/storage/homefs/jb25u082/logs/plyCheck_Project1000.err
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00:00
#SBATCH --partition epyc2
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=judit.burgayaventura@unibe.ch

echo "Launched at $(date)"
echo "Job ID: ${SLURM_JOBID}"
echo "Node list: ${SLURM_NODELIST}"
echo "Submit dir.: ${SLURM_SUBMIT_DIR}"
echo "Numb. of cores: ${SLURM_CPUS_PER_TASK}"

module load Anaconda3
eval "$(conda shell.bash hook)"

conda activate snakemake

./snakemake.sh

conda deactivate
