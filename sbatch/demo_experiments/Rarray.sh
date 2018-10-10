#!/bin/bash

module load R/3.5.1-mkl

sbatch --export arg1=2 --array=0-10%3 /home/kevinl1/selectivemodel/sbatch/demo_experiments/Rarray.sbatch
