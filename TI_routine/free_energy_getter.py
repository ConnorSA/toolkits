import lammps_python_interface as lpi
import os
import subprocess
import numpy as np
import ase,ase.io

from ase.units import kB


def write_params(TIP : lpi.ThermoIntParams, TS):
    with open(f'{TIP.toploc}/params.txt' , "w") as f:
        f.write('Paramters used: \n' +\
        f'pressure: {TIP.pressure}  \n ' +\
        f'TS temperature range: {TS.min()} - {TS.max()}  \n' +\
        f'TS steps: {len(TS)}  \n' +\
        f'number of steps: {TIP.nstep}  \n' +\
        f'thermostat: {TIP.thermostat}  \n' +\
        f'barostat: {TIP.barostat}  \n' +\
        f'averaging: {TIP.averaging}  \n' +\
        f'averaging setup: {TIP.averaging_setup}  \n' +\
        f'Number of atoms: {TIP.N}  \n' +\
        f'TS array: {TS}  \n'
        )
    return


HBAR=6.582119569e-16 ## eV per Hz
EV_TO_J=1.60218e-19
AA_TO_M=10e-10
AMU_TO_KG=1.66054e-27
BAR_TO_EV_PER_AA3=1/(160.21766208*10000)

def get_EC_FFE(k, m, T, N):
    k=k*EV_TO_J/(AA_TO_M**2) #eV/ \AA**2 to J
    m=m*AMU_TO_KG # g per mole to kg
    omega=np.sqrt(k/m) # Hz
    return 3*N*kB*T*np.log(omega*HBAR/(kB*T)) #eV


def get_COM_correction(k, N, V, T):
    return kB*T*np.log((N/V)*(2*np.pi*kB*T/(N*k))**(3/2)  )


# def get_W_irr(file, timestep, thermoprint):
#     RL=lpi.ReadLammps(file,block=1)
#     dlambda_ds = np.array(RL.dict['f_adiabat[2]'], dtype=float)
#     dlambda_dt = dlambda_ds/timestep
#     Hamiltonian = np.array(RL.dict['PotEng'], dtype=float) ##eV
#     Lambda=np.array(RL.dict['f_adiabat[1]'], dtype=float) 
#     dlambda_dt = (dlambda_dt[1:] + dlambda_dt[0:-1]) /2
#     dHamiltonian = Hamiltonian[1:] - Hamiltonian[0:-1]
#     dlambda=Lambda[1:] - Lambda[0:-1]
#     tosum = dlambda_dt*dHamiltonian/dlambda
#     W_irr_forward= np.ma.masked_invalid(tosum).sum()*timestep
#     #W_irr_forward = sum(dlambda_dt*dHamiltonian/dlambda)*timestep
#     RL=lpi.ReadLammps(file,block=3)
#     dlambda_ds = np.array(RL.dict['f_adiabat[2]'], dtype=float)
#     dlambda_dt = dlambda_ds/timestep
#     Hamiltonian = np.array(RL.dict['PotEng'], dtype=float) ##eV
#     Lambda=np.array(RL.dict['f_adiabat[1]'], dtype=float) 
#     dlambda_dt = (dlambda_dt[1:] + dlambda_dt[0:-1]) /2
#     dHamiltonian = Hamiltonian[1:] - Hamiltonian[0:-1]
#     dlambda=Lambda[1:] - Lambda[0:-1]
#     tosum = dlambda_dt*dHamiltonian/dlambda
#     W_irr_backward= np.ma.masked_invalid(tosum).sum()*timestep
#     #W_irr_backward = sum(dlambda_dt*dHamiltonian/dlambda)*timestep
#     W_irr = 0.5*(W_irr_forward - W_irr_backward)
#     return W_irr*thermoprint


def get_W_irr3(file, thermoprint):
    RL=lpi.ReadLammps(file,block=1)
    Hamiltonian_potential = np.array(RL.dict['PotEng'], dtype=float) ##eV
    Hamiltonian_springs = -np.array(RL.dict['f_adiabat'], dtype=float) #sign correction for potential energy
    Hamiltonian=Hamiltonian_potential+Hamiltonian_springs
    Lambda=np.array(RL.dict['f_adiabat[1]'], dtype=float)
    fw = np.trapz(Hamiltonian, Lambda)
    RL=lpi.ReadLammps(file,block=3)
    Hamiltonian_potential = np.array(RL.dict['PotEng'], dtype=float) ##eV
    Hamiltonian_springs = -np.array(RL.dict['f_adiabat'], dtype=float) #sign correction for potential energy
    Hamiltonian=Hamiltonian_potential+Hamiltonian_springs
    Lambda=np.array(RL.dict['f_adiabat[1]'], dtype=float)
    bw = np.trapz(Hamiltonian, Lambda)
    w = 0.5*(fw - bw)
    return w



def get_free_energy(TIP: lpi.ThermoIntParams, logfilename='thermoint.log'):
    ###set up and run setup.
    TIP.setup_crystal()
    subprocess.run(f"{TIP.LAMMPS_RUN_COMMAND} < {TIP.location}/in.setup_crystal > {TIP.location}/setup_crystal.out", shell=True)
    setup_crystal = lpi.ReadLammps(f'{TIP.location}/setup_crystal.out')
    ###average dimensions
    atoms_crystal_end = ase.io.read(f'{TIP.location}/data.setup_crystal_end', format='lammps-data', style='atomic')
    Lx, Ly, Lz = setup_crystal.average_lx(TIP.averaging_setup), setup_crystal.average_ly(TIP.averaging_setup), setup_crystal.average_lz(TIP.averaging_setup)
    atomsforTI=TIP.atoms_unrattled #defaults to equilibrated otherwise.
    atomsforTI.set_cell([Lx, Ly, Lz], scale_atoms=True)
    ase.io.write(f'{TIP.location}/data.setup_crystal_end_meaned', atomsforTI, format='lammps-data')
    # TIP.setup_getspring()
    # subprocess.run(f"{TIP.LAMMPS_RUN_COMMAND} < {TIP.location}/in.get_spring > {TIP.location}/get_spring.out", shell=True)
    # getspring_crystal = lpi.ReadLammps(f'{TIP.location}/get_spring.out')
    ### setup Einstein crystal and run Thermo Int.
    # msd = getspring_crystal.average_msd() ## \AA**2
    msd = setup_crystal.average_msd(averaging=TIP.averaging_setup) ## \AA**2
    spring_const= 3*kB*TIP.temperature/msd ## eV/ \AA**2
    TIP.setup_thermoint(spring_const)
    subprocess.run(f"{TIP.LAMMPS_RUN_COMMAND} < {TIP.location}/in.thermoint > {TIP.location}/thermoint.out", shell=True)
    thermoint_file=f"{TIP.location}/thermoint.out"
    W_irr = get_W_irr3(thermoint_file, TIP.thermoprint_switch)
    EC_ffe = get_EC_FFE(k=spring_const, m=TIP.mass, T=TIP.temperature, N=TIP.N)
    RL=lpi.ReadLammps(thermoint_file,block=0) #read first NVT of run.
    COM_correction = get_COM_correction(spring_const, TIP.N, RL.average_vol(TIP.averaging_setup), TIP.temperature)
    F0 = EC_ffe+ W_irr + COM_correction
    PV = BAR_TO_EV_PER_AA3*TIP.pressure*RL.average_vol(TIP.averaging)
    G0 = F0 + PV
    if not(os.path.exists(logfilename)):
        with open(logfilename , "w") as f:
            f.write(r'Pressure(NPT)  Pressure(NVT)  Temperature  G0/atom     F0/atom     W_irr/atom   spring   COM_correct/atom' +'\n')
    with open(logfilename , "a") as f:
        f.write(f"{TIP.pressure:13f}  {RL.average_pressure(TIP.averaging):13f}  {TIP.temperature:11f}  {G0/TIP.N:10f}  {F0/TIP.N:10f}  {W_irr/TIP.N:10f}   {spring_const:6f}   {COM_correction/TIP.N:14f}\n")
    return G0, F0, W_irr



def get_RS_path(TIP: lpi.ThermoIntParams, finaltemp):
    ###set up and run setup.
    #TIP.setup_crystal()
    #subprocess.run(f"{TIP.LAMMPS_RUN_COMMAND} < {TIP.location}/in.setup_crystal > {TIP.location}/setup_crystal.out", shell=True)
    setup_crystal = lpi.ReadLammps(f'{TIP.location}/setup_crystal.out')
    ###average dimensions
    #atoms_crystal_end = ase.io.read(f'{TIP.location}/data.setup_crystal_end', format='lammps-data', style='atomic')
    #Lx, Ly, Lz = setup_crystal.average_lx(TIP.averaging_setup), setup_crystal.average_ly(TIP.averaging_setup), setup_crystal.average_lz(TIP.averaging_setup)
    #atomsforTI=TIP.atoms_unrattled #defaults to equilibrated otherwise.
    #atomsforTI.set_cell([Lx, Ly, Lz], scale_atoms=True)
    #ase.io.write(f'{TIP.location}/data.setup_crystal_end_meaned', atomsforTI, format='lammps-data')
    ### setup Einstein crystal and run Thermo Int.
    msd = setup_crystal.average_msd(averaging=TIP.averaging_setup) ## \AA**2
    spring_const= 3*kB*TIP.temperature/msd ## eV/ \AA**2
    #spring_const=1
    lambdaf=TIP.temperature/finaltemp
    TIP.setup_RS(spring=spring_const, lambdaf=lambdaf)
    subprocess.run(f"{TIP.LAMMPS_RUN_COMMAND} < {TIP.location}/in.RS > {TIP.location}/RS.out", shell=True)
    return



