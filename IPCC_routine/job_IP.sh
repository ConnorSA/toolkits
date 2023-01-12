#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --time=12:00:00

module purge
module load iccifort/2020.4.304
module load impi/2019.9.304
module load imkl/2020.4.304
module load ASE/3.21.1


ulimit -s unlimited

python phase_diagram_via_isobars.py



