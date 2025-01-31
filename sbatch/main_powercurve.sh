module load R/3.5.1-mkl

#arguments:
## - 1st: type of signal (1 = 1-jump, 2 = 2-jump)
## - 2nd: signal to noise ratio, integers from 1 to 6 for (0, 0.25, 0.5, 1, 2, 4)
## - 3th: method (1 = bs, 2 = fl)
## - 4th: sigma (either 1 or 2 or 1 or NA)
## - 5th: ksteps (2 to 4)
## - 6th: decluttered (0 = no, 1 = yes)
## - array: index of hypothesis test

sbatch --export arg1=2,arg2=1,arg3=1,arg4=2,arg5=2,arg6=0 --array=1-500%10 /home/kevinl1/selectivemodel/sbatch/main_powercurve.sbatch
sbatch --export arg1=2,arg2=2,arg3=1,arg4=2,arg5=2,arg6=0 --array=1-500%10 /home/kevinl1/selectivemodel/sbatch/main_powercurve.sbatch
sbatch --export arg1=2,arg2=3,arg3=1,arg4=2,arg5=2,arg6=0 --array=1-500%10 /home/kevinl1/selectivemodel/sbatch/main_powercurve.sbatch
sbatch --export arg1=2,arg2=4,arg3=1,arg4=2,arg5=2,arg6=0 --array=1-500%10 /home/kevinl1/selectivemodel/sbatch/main_powercurve.sbatch
sbatch --export arg1=2,arg2=5,arg3=1,arg4=2,arg5=2,arg6=0 --array=1-500%10 /home/kevinl1/selectivemodel/sbatch/main_powercurve.sbatch
sbatch --export arg1=2,arg2=6,arg3=1,arg4=2,arg5=2,arg6=0 --array=1-500%10 /home/kevinl1/selectivemodel/sbatch/main_powercurve.sbatch
