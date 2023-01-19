from lammps_input_writer_EAM import *
from calculate_mu import *
import ase, ase.io, ase.build
import os
from warnings import simplefilter
simplefilter("ignore")


#
mpirun='/warwick/desktop/2018/software/Compiler/intel/2017.4.196-GCC-6.4.0-2.28/impi/2017.3.196/bin64/mpirun'
# e.g. "mpirun -np 4 ./lmp_mpi" or "srun ./lmp_mpiicpc"
LAMMPS_RUN_COMMAND = f"{mpirun} -np 20 ./lmp_mpi"


#run folder
folder=f'test'
if not(os.path.exists(folder)):
    os.mkdir(folder)


#building things
topdir=os.getcwd()
atoms=ase.io.read('bcc_Ti_geomopti.castep')
a = atoms.cell.cellpar()[0]*np.sqrt(4/3)
atoms = ase.build.bulk('Ti', crystalstructure='bcc', a=a, orthorhombic='True')
N_x=6
N_y=6
N_z=20
sc=[N_x,N_y,N_z]
atoms = atoms*sc
atoms.rattle(0.10)


IP=IPparams(toploc=folder,
            pressure=10000,
            temperature=1500,
            step=25000,
            hkl=[N_x,N_y,0],
            thermostat=0.25,
            barostat=1.00,
	        Nz=N_z,
            thinned=30,
            atoms=atoms.copy(),
            LAMMPS_RUN_COMMAND=LAMMPS_RUN_COMMAND,
            auto=True)
IP.write_IP_params()


#run one IP routine.
run_till_converged(IP, IP.traj_name, melt_steps=10, tol=1, samples=5)


#phase diagram via isobars
pressures=[IP.pressure,20000,30000]
for i, p in enumerate(pressures[1:]):
    crystal_cc = ReadLammps(f'{IP.location}/crystal_auto_EAM.out')
    liquid_cc = ReadLammps(f'{IP.location}/liquid_auto_EAM.out')
    next_temp=classius_clapeyron_next_temp(current_T=IP.temperature, dP=p-pressures[i], crystal=crystal_cc, liquid=liquid_cc, thinned=IP.thinned)
    IP.next_isobar_start(pressure=p, temperature=next_temp)
    run_till_converged(IP, IP.traj_name, melt_steps=10, tol=1, samples=5)


