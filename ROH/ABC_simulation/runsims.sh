#!/usr/bin/bash
#SBATCH -p long
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-63
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load R
module load clusterbasics

php runsim.php ${SLURM_ARRAY_TASK_COUNT}  ${SLURM_ARRAY_TASK_ID}  1000000 
