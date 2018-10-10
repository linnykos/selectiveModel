#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=28
#SBATCH -t 8:00:00
#SBATCH --mem 40000
#SBATCH --job-name=main_powercurve
#SBATCH --output=slurm_%a.out

/opt/packages/R/3.5.1-mkl/lib64/R/bin/Rscript --vanilla /home/kevinl1/selectivemodel/sbatch/main_powercurve.R $arg1 $arg2 $arg3 $arg4 $arg5 $arg6 $arg7
