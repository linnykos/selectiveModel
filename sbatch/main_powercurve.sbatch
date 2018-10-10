#!/bin/bash
#SBATCH --job-name=main_powercurve
#SBATCH --output=/home/kevinl1/selectivemodel/sbatch/logs/slurm.%A.%a.%j.out
#SBATCH --error=/home/kevinl1/selectivemodel/sbatch/logs/slurm.%A.%a.%j.out
#SBATCH --time=00:10:00
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks=1

/opt/packages/R/3.5.1-mkl/lib64/R/bin/Rscript --vanilla /home/kevinl1/selectivemodel/sbatch/main_powercurve.R $arg1 $arg2 $arg3 $arg4 $arg5 $arg6 ${SLURM_ARRAY_TASK_ID}
