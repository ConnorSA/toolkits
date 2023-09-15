import ase, ase.build, ase.io
import numpy
from lammps_python_interface import ThermoIntParams
from free_energy_getter import get_free_energy, get_RS_path
import os
import numpy as np
#LAMMPS_RUN_COMMAND = "mpirun -np 16 ./lmp_g++_openmpi"
LAMMPS_RUN_COMMAND = f"srun ./lmp_mpiicpc"


atoms=ase.build.bulk('Ti', orthorhombic=True)#, crystalstructure='bcc', a=3.2)
sc=[10,10,10]
atoms=atoms*sc
atoms_unrattled=atoms.copy()
atoms.rattle(0.05)

Trange=np.linspace(800,1300,51)
#potential_str=("pair_style      quip \n"
#               "pair_coeff      * * Ti_turbo_awfr_0.02_min_0.025_cwvr_0.05_min_0.10_SP3000_curpoints.xml \"Potential xml_label=GAP_2023_8_10_60_16_29_2_871\" 22 \n")

pair_style = "eam/fs"
pair_coeff = "* * Ti1.eam.fs Ti"

folder=f'RS-HCP-{sc[0]}x{sc[1]}x{sc[2]}-{Trange[0]}-{Trange[-1]}-longboi-swapped'
if not(os.path.exists(folder)):
    os.mkdir(folder)

TIP=ThermoIntParams(
    LAMMPS_RUN_COMMAND=LAMMPS_RUN_COMMAND, 
    atoms=atoms,
    atoms_unrattled=atoms_unrattled, 
    toploc=folder, 
    pair_style=pair_style,
    pair_coeff=pair_coeff,
    mass=48,
    pressure=0, 
    temperature=Trange[0], 
    timestep=0.002, 
    nstep=200000,
    nstep_setup=50000,
    nstep_eq=50000,  
    thermostat=0.95, 
    barostat=2.00,
    averaging_setup=100,
    averaging=1000,
    thermoprint=50
)


get_free_energy(TIP, logfilename=f'thermoint_{folder}.log')
get_RS_path(TIP, finaltemp=1300)

#for T in Trange:
#    TIP.update_PT(new_pressure=0, new_temperature=T)
#    get_free_energy(TIP, logfilename=f'thermoint_{folder}.log')
