#!/bin/bash
#SBATCH --job-name=RarrayJob
#SBATCH --output=/home/kevinl1/selectivemodel/sbatch/logs/slurm.%A.%a.%j.out
#SBATCH --error=/home/kevinl1/selectivemodel/sbatch/logs/slurm.%A.%a.%j.out
#SBATCH --time=00:05:00
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks=1

/opt/packages/R/3.5.1-mkl/lib64/R/bin/Rscript --vanilla nodename.R $arg1 ${SLURM_ARRAY_TASK_ID}
