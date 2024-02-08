#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --time=24:00:00

module purge
module load iccifort/2020.4.304 impi/2019.9.304 imkl/2020.4.304 ASE/3.21.1


ulimit -s unlimited
JOB_NAME=${SLURM_JOB_NAME}


python driver_RS.py ${JOB_NAME}



