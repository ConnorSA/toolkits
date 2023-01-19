#!/bin/bash

module purge
module load GCC/8.3.0 OpenBLAS/0.3.7

MPI_COMMAND=/warwick/desktop/2018/software/Compiler/intel/2017.4.196-GCC-6.4.0-2.28/impi/2017.3.196/bin64/mpirun

python3 phase_diagram_via_isobars_clean.py
