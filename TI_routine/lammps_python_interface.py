import numpy as np
import io, os
import ase, ase.io
from ase.units import kB


### Notes to self
# Probably want to do the NPT side of simulations consistently. Authors used NPH + Langevin, and constrained COM.
# I think this means that when they do the thermodynamic integration, it's not actually in NVT?

###    f"fix ensemble all npt temp {temperature} {temperature} {thermostat} {isostring} {pressure} {pressure} {barostat} pchain 32  \n" \

    # f"variable xcm equal xcm(all,x) \n"\
    # f"variable ycm equal xcm(all,y) \n"\
    # f"variable zcm equal xcm(all,z) \n"\
    # f"fix f1 all nph aniso  {pressure} {pressure}  {barostat} fixedpoint" +r" ${xcm} ${ycm} ${zcm}" + "\n"\
    # f"fix f2 all langevin {temperature} {temperature} {thermostat} 999 zero yes \n" \
    # f"compute c1 all temp/com \n" \
    # f"fix_modify f1 temp c1 \n" \
    # f"fix_modify f2 temp c1 \n" \

    # f"fix ensemble all npt temp {temperature} {temperature} {thermostat} {isostring} {pressure} {pressure} {barostat} pchain 32  \n" \

def write_setup_crytal(pair_style : str, pair_coeff : str, location: str, mass : float, pressure : float, temperature : float, timestep : float ,nstep = 10000,
    barostat=1.00, thermostat=0.25, aniso=True, thermoprint=50):
    isostring='aniso'
    if aniso==False:
        isostring='iso'
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.crystal  \n" \
    f"mass          1 {mass}  \n" \
    f"pair_style  {pair_style} \n" \
    f"pair_coeff  {pair_coeff} \n" \
    "neighbor        0.3 bin  \n" \
    "neigh_modify    delay 10  \n" \
    f"timestep        {timestep}  \n" \
    "run_style	verlet  \n" \
    f"velocity all create  {temperature} 1 mom yes rot yes dist gaussian \n" \
    f"variable xcm equal xcm(all,x) \n"\
    f"variable ycm equal xcm(all,y) \n"\
    f"variable zcm equal xcm(all,z) \n"\
    f"fix f1 all nph aniso  {pressure} {pressure}  {barostat} fixedpoint" +r" ${xcm} ${ycm} ${zcm}" + "\n"\
    f"fix f2 all langevin {temperature} {temperature} {thermostat} 999 zero yes \n" \
    f"compute c1 all temp/com \n" \
    f"fix_modify f1 temp c1 \n" \
    f"fix_modify f2 temp c1 \n" \
    f"compute 1 all msd com yes average yes \n" \
    f"variable msd equal c_1[4] \n" \
    f"thermo {thermoprint}  \n" \
    "thermo_style custom step temp press c_1[4] vol etotal enthalpy lx ly lz density \n" \
    f"dump dumpXYZ all xyz 100 {location}/traj_crystal.xyz  \n" \
    f"run {nstep}  \n" \
    f"write_data {location}/data.setup_crystal_end  \n"
    with open(f'{location}/in.setup_crystal' , "w") as f:
        f.write(param_string)
    return


def write_getspring(pair_style : str, pair_coeff : str, location: str, mass : float, temperature : float, timestep : float ,nstep = 10000, nstep_eq=10000,
    thermostat=0.25, aniso=True, thermoprint=10):
    isostring='aniso'
    if aniso==False:
        isostring='iso'
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.setup_crystal_end_meaned  \n" \
    f"mass          1 {mass}  \n" \
    f"pair_style  {pair_style} \n" \
    f"pair_coeff  {pair_coeff} \n" \
    "#neighbor        0.3 bin  \n" \
    "#neigh_modify    delay 10  \n" \
    f"timestep        {timestep}  \n" \
    "run_style	verlet  \n" \
    f"variable xcm equal xcm(all,x) \n"\
    f"variable ycm equal xcm(all,y) \n"\
    f"variable zcm equal xcm(all,z) \n"\
    f"fix f1 all nvt temp {temperature} {temperature} {thermostat} fixedpoint" +r" ${xcm} ${ycm} ${zcm}" + "\n"\
    f"compute  Tcm all temp/com \n" \
    f"fix_modify        f1 temp Tcm \n"\
    f"velocity all create  {temperature} 1 mom yes rot yes dist gaussian \n" \
    f"compute 1 all msd com yes average yes \n" \
    f"variable msd equal c_1[4] \n" \
    f"thermo {thermoprint}  \n" \
    "thermo_style custom step temp press c_1[4] vol etotal density \n" \
    f"dump dumpXYZ all xyz 100 {location}/traj_getspring.xyz  \n" \
    f"run {nstep} \n"
    with open(f'{location}/in.get_spring' , "w") as f:
        f.write(param_string) 
    return


### COM adjust NPT
    # f"velocity all create  {temperature} 1 mom yes rot yes \n" \
    # f"variable xcm equal xcm(all,x) \n"\
    # f"variable ycm equal xcm(all,y) \n"\
    # f"variable zcm equal xcm(all,z) \n"\
    # f"fix f1 all nph aniso  {pressure} {pressure}  {barostat} fixedpoint" +r" ${xcm} ${ycm} ${zcm}" + "\n"\
    # f"fix f2 all langevin {temperature} {temperature} {thermostat} 999 zero yes \n" \
    # f"compute c1 all temp/com \n" \
    # f"fix_modify f1 temp c1 \n" \
    # f"fix_modify f2 temp c1 \n" \

    # f"variable xcm equal xcm(all,x) \n"\
    # f"variable ycm equal xcm(all,y) \n"\
    # f"variable zcm equal xcm(all,z) \n"\
#    f"fix f1 all nvt temp {temperature} {temperature} {thermostat} fixedpoint" +r" ${xcm} ${ycm} ${zcm}" + "\n"\
    # f"velocity all create  {temperature} 1 mom yes rot yes \n" \
    # f"compute c1 all temp/com \n" \
    # f"fix_modify f1 temp c1 \n" \
def write_thermoint(pair_style : str, pair_coeff : str, location: str, mass : float, temperature : float, pressure : float, timestep : float, spring ,nstep = 10000, nstep_eq=10000,
    barostat=1.00, thermostat=0.25, aniso=True, thermoprint=50):
    isostring='aniso'
    if aniso==False:
        isostring='iso'
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.setup_crystal_end_meaned  \n" \
    f"mass          1 {mass}  \n" \
    f"pair_style  {pair_style} \n" \
    f"pair_coeff  {pair_coeff} \n" \
    "#neighbor        0.3 bin  \n" \
    "#neigh_modify    delay 10  \n" \
    f"timestep        {timestep}  \n" \
    "run_style	verlet  \n" \
    f"fix f0 all nve \n" \
    f"fix adiabat all ti/spring {spring} {nstep} {nstep_eq} function 2 \n" \
    f"fix f1 all langevin {temperature} {temperature} {thermostat} 999 zero yes \n" \
    f"compute  Tcm all temp/com \n" \
    f"fix_modify        f1 temp Tcm \n"\
    f"velocity all create  {temperature} 1 mom yes rot yes dist gaussian \n" \
    f"thermo {thermoprint}  \n" \
    f"thermo_style custom step temp press f_adiabat f_adiabat[1] f_adiabat[2] vol etotal pe enthalpy lx ly lz density \n"\
    f"dump dumpXYZ all xyz 100 {location}/traj_thermoint.xyz  \n" \
    f"run {nstep_eq} \n" \
    f"run {nstep} \n" \
    f"run {nstep_eq} \n" \
    f"run {nstep} \n" 
    with open(f'{location}/in.thermoint' , "w") as f:
        f.write(param_string) 
    return

def write_thermoint2(pair_style : str, pair_coeff : str, location: str, temperature : float, pressure : float,mass : float, timestep : float, spring ,nstep = 10000, nstep_eq=10000,
    barostat=1.00, thermostat=0.25, aniso=True, thermoprint=50):
    isostring='aniso'
    if aniso==False:
        isostring='iso'
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.setup_crystal_end_meaned  \n" \
    f"mass          1 {mass}  \n" \
    f"pair_style  {pair_style} \n" \
    f"pair_coeff  {pair_coeff} \n" \
    "neighbor        0.3 bin  \n" \
    "neigh_modify    delay 10  \n" \
    f"timestep        {timestep}  \n" \
    "run_style	verlet  \n" \
    f"fix adiabat all ti/spring {spring} {nstep} {nstep_eq} \n" \
    f"fix f1 all nvt temp {temperature} {temperature} {thermostat} \n"\
    f"velocity all create  {temperature} 1 mom yes rot yes dist gaussian \n" \
    f"compute c1 all temp/com \n" \
    f"thermo {thermoprint}  \n" \
    f"thermo_style custom step temp press f_adiabat f_adiabat[1] f_adiabat[2] vol etotal pe enthalpy lx ly lz density c_c1 \n"\
    r"#thermo_modify format f_adiabat[2] %.8f "+ "\n"\
    f"dump dumpXYZ all xyz 100 {location}/traj_thermoint.xyz  \n" \
    f"run {nstep_eq} \n" \
    f"run {nstep} \n" \
    f"run {nstep_eq} \n" \
    f"run {nstep} \n" 
    with open(f'{location}/in.thermoint' , "w") as f:
        f.write(param_string) 
    return

### Maybe don't need spring at all?
#    f"fix adiabat all ti/spring {spring} {nstep} {nstep_eq} \n" \


#    f"fix f1 all nph aniso  {pressure} {pressure}  {barostat} fixedpoint" +r" ${xcm} ${ycm} ${zcm}" + "\n"\
#    f"fix f2 all langevin {temperature} {temperature} {thermostat} 999 zero yes \n" \


def write_RS(pair_style : str, pair_coeff : str, location: str, mass : float, temperature : float, pressure : float,
              timestep : float, spring, lambdaf, nstep = 10000, nstep_eq=10000,
              barostat=1.00, thermostat=0.25, aniso=True):
    isostring='aniso'
    if aniso==False:
        isostring='iso'
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.setup_crystal_end_meaned  \n" \
    f"mass          1 {mass}  \n" \
    f"pair_style  {pair_style} \n" \
    f"pair_coeff  {pair_coeff} \n" \
    "neighbor        0.3 bin  \n" \
    "neigh_modify    delay 10  \n" \
    f"timestep        {timestep}  \n" \
    "run_style	verlet  \n" \
    f"variable lambda equal 1 \n"\
    f"fix f3 all adapt 1 pair {pair_style} scale * * v_lambda  \n"\
    f"variable xcm equal xcm(all,x) \n"\
    f"variable ycm equal xcm(all,y) \n"\
    f"variable zcm equal xcm(all,z) \n"\
    f"fix f2 all langevin {temperature} {temperature} {thermostat} 999 zero yes \n"\
    f"fix f1 all nph aniso  {pressure} {pressure}  {barostat} fixedpoint" +r" ${xcm} ${ycm} ${zcm}" + "\n"\
    f"compute c1 all temp/com \n" \
    f"fix_modify f1 temp c1 \n" \
    f"fix_modify f2 temp c1 \n" \
    f"velocity all create  {temperature} 1 mom yes rot yes \n" \
    f"thermo 10  \n" \
    f"thermo_style custom step temp press v_lambda vol etotal pe enthalpy lx ly lz density c_c1 \n"\
    r"#thermo_modify format f_adiabat[2] %.8f "+ "\n"\
    f"dump dumpXYZ all xyz 100 {location}/traj_RS.xyz  \n" \
    f"run {nstep_eq} \n" \
    f"variable lambda equal 1/(1+elapsed/{nstep}*(1/{lambdaf}-1)) \n"\
    f"run {nstep} \n" \
    f"variable lambda equal {lambdaf} \n"\
    f"run {nstep_eq} \n" \
    f"variable lambda equal 1/(1+(1-elapsed/{nstep})*(1/{lambdaf}-1)) \n"\
    f"run {nstep} \n" 
    with open(f'{location}/in.RS' , "w") as f:
        f.write(param_string) 
    return



class ThermoIntParams():
    def __init__(self, LAMMPS_RUN_COMMAND, atoms, toploc, pair_style, pair_coeff ,pressure, temperature, mass, timestep, nstep,  
    thermostat, barostat, nstep_setup=None,nstep_eq=None, averaging_setup=30,averaging=30, aniso=True, atoms_unrattled=None, thermoprint=50):
        self.toploc = toploc
        self.location=f'{toploc}/P_{str(pressure).zfill(6)}_T_{str(temperature).zfill(4)}'
        if not(os.path.exists(self.location)):
            os.mkdir(self.location)
        ase.io.write(f'{self.location}/data.crystal', atoms, format='lammps-data')
        self.pair_style=pair_style
        self.pair_coeff=pair_coeff
        self.pressure=pressure
        self.mass = mass
        self.barostat=barostat
        self.temperature=temperature
        self.thermostat=thermostat
        self.timestep=timestep
        self.nstep=nstep
        self.aniso=aniso
        if nstep_eq==None:
            self.nstep_eq=nstep
        else:    
            self.nstep_eq=nstep_eq

        if nstep_setup==None:
            self.nstep_setup=nstep
        else:    
            self.nstep_setup=nstep_setup
        self.averaging_setup=averaging_setup
        self.averaging=averaging
        self.N = len(atoms)
        self.atoms = atoms
        if not(atoms_unrattled ==None):
            self.atoms_unrattled = atoms_unrattled
        else: self.atoms_unrattled = atoms
        self.LAMMPS_RUN_COMMAND = LAMMPS_RUN_COMMAND
        self.traj_name =  f"{toploc}/P_{str(pressure).zfill(6)}.iptraj"
        self.thermoprint=thermoprint
        self.thermoprint_switch=10 #this needs to be very small for converged integration
        return


    def update_PT(self, new_pressure, new_temperature):
        self.location=f'{self.toploc}/P_{str(new_pressure).zfill(6)}_T_{str(new_temperature).zfill(4)}'
        if not(os.path.exists(self.location)):
            os.mkdir(self.location)
        ase.io.write(f'{self.location}/data.crystal', self.atoms, format='lammps-data')
        self.pressure=new_pressure
        self.temperature=new_temperature
        return

    def setup_crystal(self):
        write_setup_crytal(pair_style=self.pair_style, pair_coeff=self.pair_coeff, location=self.location, mass=self.mass,pressure=self.pressure,
                           temperature=self.temperature, 
                           timestep=self.timestep, nstep=self.nstep_setup,
                           barostat=self.barostat, thermostat=self.thermostat, aniso=self.aniso, thermoprint=self.thermoprint)
        return

    def setup_getspring(self):
        write_getspring(pair_style=self.pair_style, pair_coeff=self.pair_coeff, location=self.location, mass=self.mass,
                           temperature=self.temperature, 
                           timestep=self.timestep, nstep=self.nstep_eq, thermostat=self.thermostat, aniso=self.aniso,  thermoprint=self.thermoprint)
        return


    def setup_thermoint(self, spring):
        write_thermoint(pair_style=self.pair_style, pair_coeff=self.pair_coeff, location=self.location, mass=self.mass, temperature=self.temperature,
                          pressure=self.pressure, timestep=self.timestep,spring=spring , nstep=self.nstep, nstep_eq=self.nstep_eq,
                           barostat=self.barostat, thermostat=self.thermostat, aniso=self.aniso, thermoprint=self.thermoprint_switch)
        return
    
    def setup_RS(self, spring, lambdaf):
        write_RS(pair_style=self.pair_style, pair_coeff=self.pair_coeff, location=self.location, mass=self.mass, temperature=self.temperature,
                          pressure=self.pressure, timestep=self.timestep,spring=spring,lambdaf=lambdaf , nstep=self.nstep, nstep_eq=self.nstep_eq,
                           barostat=self.barostat, thermostat=self.thermostat, aniso=self.aniso)
        return

class ReadLammps():
    #Read lammps output: assumes Step is outputted.
    def __init__(self, filein,block=0):
        self.block=block
        self.get_dict(filein)
        
        
    def get_dict(self,filein):
        b=0
        readcheck=False
        with open(filein) as f:
            lines=f.readlines()
        self.lines=lines
        for i, val in enumerate(lines):
            if 'Step' in val and b==self.block:
                readcheck=True
                break
            elif 'Step' in val:
                b+=1
        if readcheck == False: print('EOF: did not parse properly')

        start_loc=i+1      
        for i, val in enumerate(lines[start_loc:]):
            if 'Loop time' in val:
                break
        end_loc=i+start_loc
        labels = lines[start_loc-1].split()[0:]
        mydict = {}
        for i in labels:
            mydict[i]=[]
        for i, val in enumerate(lines[start_loc:end_loc]):
                for j, val2 in enumerate(labels):
                    mydict[val2].append(val.split()[j])
        self.dict = mydict

    def average_vol(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Volume'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Volume'], dtype=float).mean()

    def average_density(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Density'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Density'], dtype=float).mean()

    def average_energy(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['TotEng'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['TotEng'], dtype=float).mean()


    def average_temp(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Temp'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Temp'], dtype=float).mean()

    def average_pressure(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Press'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Press'], dtype=float).mean()

    def std_pressure(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Press'], dtype=float)[-averaging:].std()
        return np.array(self.dict['Press'], dtype=float).std()

    def average_enthalpy(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Enthalpy'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Enthalpy'], dtype=float).mean()

    def average_lx(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Lx'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Lx'], dtype=float).mean()

    def average_ly(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Ly'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Ly'], dtype=float).mean()

    def average_lz(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['Lz'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['Lz'], dtype=float).mean()
    
    def average_msd(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['c_1[4]'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['c_1[4]'], dtype=float).mean()
    
    def average_EC_energy(self, averaging=None):
        if averaging!=None and type(averaging) == int:
            return np.array(self.dict['f_adiabat'], dtype=float)[-averaging:].mean()
        return np.array(self.dict['f_adiabat'], dtype=float).mean()
