module load R/3.5.1-mkl

#arguments:
## - 1st: type of signal (1 = 1-jump, 2 = 2-jump)
## - 2nd: signal to noise ratio, integers from 1 to 6 for (0, 0.25, 0.5, 1, 2, 4)
## - 3th: method (1 = bs, 2 = fl)
## - 4th: sigma (either 1 or 2 or 1 or NA)
## - 5th: ksteps (2 to 4)
## - 6th: decluttered (0 = no, 1 = yes)
## - 7th: trials

sbatch --export arg1=2,arg2=6,arg3=1,arg4=2,arg5=2,arg6=0,arg7=5,arg8=1 -o /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out -e /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out /home/kevinl1/selectivemodel/sbatch/demo_experiments/main_powercurve_1_for.sbatch
sbatch --export arg1=2,arg2=6,arg3=1,arg4=2,arg5=2,arg6=0,arg7=5,arg8=1 -o /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out -e /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out /home/kevinl1/selectivemodel/sbatch/demo_experiments/main_powercurve_1_mlapply1.sbatch
sbatch --export arg1=2,arg2=6,arg3=1,arg4=2,arg5=2,arg6=0,arg7=5,arg8=28 -o /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out -e /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out /home/kevinl1/selectivemodel/sbatch/demo_experiments/main_powercurve_28_for.sbatch
sbatch --export arg1=2,arg2=6,arg3=1,arg4=2,arg5=2,arg6=0,arg7=5,arg8=28 -o /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out -e /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out /home/kevinl1/selectivemodel/sbatch/demo_experiments/main_powercurve_28_mlapply1.sbatch
sbatch --export arg1=2,arg2=6,arg3=1,arg4=2,arg5=2,arg6=0,arg7=5,arg8=28 -o /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out -e /home/kevinl1/selectivemodel/sbatch/logs/slurm.%N.%j.out /home/kevinl1/selectivemodel/sbatch/demo_experiments/main_powercurve_28_mlapply28.sbatch
