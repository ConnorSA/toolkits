import subprocess
from lammps_input_writer_EAM import *
from calculate_mu import *
import ase, ase.io, ase.build
import os
from warnings import simplefilter
simplefilter("ignore")




#run folder
folder=f'rundir_7x7x50_kappa_pref_1.000_very_coldstart'
if not(os.path.exists(folder)):
    os.mkdir(folder)


#building things
topdir=os.getcwd()
atoms=ase.io.read('bcc_Ti_geomopti.castep')
a = atoms.cell.cellpar()[0]*np.sqrt(4/3)
atoms = ase.build.bulk('Ti', crystalstructure='bcc', a=a, orthorhombic='True')
N_xy=7
N_z=50
sc=[N_xy,N_xy,N_z]
HKL=[N_xy,N_xy,0]
atoms = atoms*sc
atoms.rattle(0.10)
z_half=np.round(atoms.cell[2][2]/2, 1)


#IP params for manual set or auto. auto overwrites manual.
IP_a=None
IP_kappa=None
auto_IP=True
#Scale the value of kappa (for auto IP)
kappa_prefactor=1
#Scale the size of the next temp step
step_prefactor=1
#max number of IP per isobar
melt_steps=50
#number of samples applied IP output for averaging and tol check
samples=50
#tol of the IP routine in K, when to stop at curreent P.
tol=0


#thermostat and barostat
start_T = 300
start_P = 10000
barostat=1.00
thermostat=0.25
#number of MD steps
Nstep=25000
#amount of samples applied to LAMMPS output for averaging.
thinning=30





subfolder=f'/initialised_P_{str(start_P).zfill(6)}'
if not(os.path.exists(folder+subfolder)):
    os.mkdir(folder+subfolder)
ase.io.write(f'{folder+subfolder}/data.hcp_Ti_crystal', atoms, format='lammps-data')

IP=IPparams(location=folder+subfolder,
            pressure=start_P,
            temperature=start_T,
            step=Nstep,
            hkl=HKL,
            thermostat=thermostat,
            barostat=barostat,
            z_half=z_half,
	        Nz=N_z,
            a=IP_a,
            kappa=IP_kappa,
            thinned=thinning,
            N=len(atoms),
            atoms=atoms.copy())

with open(f'{folder}/params.txt' , "w") as f:
    f.write('Paramters used: \n' +\
    f'start pressure: {IP.pressure}  \n ' +\
    f'start temperature: {IP.temperature}  \n' +\
    f'number of steps: {IP.step}  \n' +\
    f'hkl: {IP.hkl}  \n' +\
    f'thermostat: {IP.thermostat}  \n' +\
    f'barostat: {IP.barostat}  \n' +\
    f'z_half: {IP.z_half}  \n' +\
    f'IP_a: {IP.a}  \n' +\
    f'IP_kappa: {IP.kappa}  \n' +\
    f'IP_kappa: {auto_IP}  \n' +\
    f'thinned: {IP.thinned}  \n' +\
    f'Number of atoms: {IP.N}  \n' +\
    f'Kappa prefactor: {kappa_prefactor}  \n' +\
    f'Next T step prefactor: {step_prefactor}  \n'
    )


subfolder=f'/melt_isobar_P_{str(start_P).zfill(6)}'
if not(os.path.exists(folder+subfolder)):
    os.mkdir(folder+subfolder)
ase.io.write(f'{folder+subfolder}/data.hcp_Ti_crystal', atoms, format='lammps-data')
IP.location = folder+subfolder

#first run
tmeltout=f'{folder}/melt_isobar_P_{str(start_P).zfill(6)}.txt'
run_till_converged(IP, tmeltout, melt_steps=melt_steps, tol=tol, samples=samples, kappa_prefactor=kappa_prefactor,
update_IP_params=auto_IP, step_prefactor=step_prefactor)

# dP = 10000 #step by 10000 bar
# for i in range(1):
#     IP.pressure = IP.pressure + dP
#     tmeltout=f'{folder}/melt_isobar_P_{str(IP.pressure).zfill(6)}.txt'
#     subfolder=f'/melt_isobar_P_{str(IP.pressure).zfill(6)}'
#     if not(os.path.exists(folder+subfolder)):
#         os.mkdir(folder+subfolder)
#     crystal_cc = ReadLammps(f'{IP.location}/crystal_auto_EAM.out')
#     liquid_cc = ReadLammps(f'{IP.location}/liquid_auto_EAM.out')
#     IP.location=folder+subfolder
#     ase.io.write(f'{folder+subfolder}/data.hcp_Ti_crystal', atoms, format='lammps-data')
#     IP.temperature=classius_clapeyron_next_temp(current_T=IP.temperature, dP=dP, crystal=crystal_cc, liquid=liquid_cc, thinned=IP.thinned)
#     run_till_converged(IP, tmeltout, melt_steps=melt_steps, tol=tol, samples=samples, kappa_prefactor=kappa_prefactor, 
#     update_IP_params=auto_IP, step_prefactor=step_prefactor)
