#!/bin/bash

#SBATCH --job-name=RarrayJob
#SBATCH --output=RarrayJob_%A_%a.out
#SBATCH --error=RarrayJob_%A_%a.err
#SBATCH --array=1-4
#SBATCH --time=00:05:00
#SBATCH --partition=RM-small
#SBATCH --nodes=1
#SBATCH --ntasks=28

export HOSTNAME=`cat /etc/hostname`

cd $SLURM_SUBMIT_DIR
mkdir -p dir-${SLURM_ARRAY_TASK_ID}
cd       dir-${SLURM_ARRAY_TASK_ID}


echo "My Hostname:            " $HOSTNAME             >  ${SLURM_ARRAY_TASK_ID}.log
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID  >> ${SLURM_ARRAY_TASK_ID}.log

Rscript --vanilla nodename.R >> ${SLURM_ARRAY_TASK_ID}.log
