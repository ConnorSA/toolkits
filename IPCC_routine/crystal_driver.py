import subprocess
from lammps_input_writer_EAM import *
from calculate_mu import *
import ase, ase.io, ase.build
import os
from warnings import simplefilter
simplefilter("ignore")


def write_in_crystal_nvt(location: str, pressure : float, temperature : float, step = 10000, hkl=[6,6,0],
    barostat=1.00, thermostat=0.25):
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.hcp_Ti_crystal  \n" \
    "mass          1 48.0  \n" \
    "pair_style eam/fs  \n" \
    "pair_coeff * * Ti1.eam.fs Ti  \n" \
    "neighbor        0.3 bin  \n" \
    "neigh_modify    delay 10  \n" \
    "timestep        0.002  \n" \
    "run_style	verlet  \n" \
    f"fix ensemble all nvt temp {temperature} {temperature} {thermostat}  \n" \
    f"fix bias all rhok {hkl[0]} {hkl[1]} {hkl[2]} 0.0 0.0  \n" \
    "thermo 50  \n" \
    "thermo_style custom step temp press density f_bias[3] vol etotal enthalpy lx ly lz\n" \
    f"run {step}  \n" 
    with open(f'{location}/in.crystal_auto_EAM' , "w") as f:
        f.write(param_string)
    return


#run folder
folder=f'rundir_7x7x50'
if not(os.path.exists(folder)):
    os.mkdir(folder)


topdir=os.getcwd()
atoms=ase.io.read('bcc_Ti_geomopti.castep')
a = atoms.cell.cellpar()[0]*np.sqrt(4/3)
atoms = ase.build.bulk('Ti', crystalstructure='bcc', a=a, orthorhombic='True')
sc=[7,7,50]
HKL=[7,7,0]
atoms = atoms*sc
atoms.rattle(0.10)


#thermostat and barostat
start=200
stop=2000
templist = np.linspace(start, stop, 1+ (stop-start)//10)
start_P = 10000
barostat=1.00
thermostat=0.25

#number of MD steps
Nstep=50000

#amount of samples applied to LAMMPS output for averaging.
thinning=50


ase.io.write(f'{folder}/data.hcp_Ti_crystal', atoms, format='lammps-data')
for temperature in templist:    
    write_in_crystal_nvt(location=folder, temperature=temperature, thermostat=thermostat, step=Nstep)    
    subprocess.run(f"srun ./lmp_mpiicpc < {folder}/in.crystal_auto_EAM > {folder}/crystal_temperature_{str(temperature).zfill(4)}.out", shell=True) 

