import numpy as np
import io
import ase, ase.io
from ase.units import kB

def write_in_crystal(location: str, pressure : float, temperature : float, step = 10000, hkl=[6,6,0],
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
    f"fix ensemble all npt temp {temperature} {temperature} {thermostat} aniso {pressure} {pressure} {barostat} pchain 32  \n" \
    f"fix bias all rhok {hkl[0]} {hkl[1]} {hkl[2]} 0.0 0.0  \n" \
    "thermo 50  \n" \
    "thermo_style custom step temp press density f_bias[3] vol etotal enthalpy lx ly lz\n" \
    f"dump dumpXYZ all xyz 500 {location}/traj_crystal_EAM.xyz  \n" \
    f"run {step}  \n" \
    f"write_data {location}/data.crystal_end_EAM  \n"
    with open(f'{location}/in.crystal_auto_EAM' , "w") as f:
        f.write(param_string)
    return 
    
def write_in_setup_pinning(location: str, pressure : float, temperature : float, step = 10000, hkl=[6,6,0], z_half=0, 
    barostat=1.00, thermostat=0.25):
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.crystal_end_meaned  \n" \
    "mass          1 48.0  \n" \
    "pair_style eam/fs  \n" \
    "pair_coeff * * Ti1.eam.fs Ti  \n" \
    "neighbor        0.3 bin  \n" \
    "timestep        0.002  \n" \
    "run_style	verlet  \n" \
    f"region left plane 0 0 {z_half} 0 0 1  \n" \
    "group left region left  \n" \
    "velocity left create 4000.0 1 mom yes rot yes  \n" \
    "fix ensemble left nve     # Note: only move particle in left-hand side  \n" \
    f"fix langevin left langevin 4000.0 {temperature} {thermostat} 2017  \n" \
    "thermo_style custom step temp pzz pe lz vol density etotal  \n" \
    "thermo 100  \n" \
    f"#dump dumpXYZ all xyz 50 {location}/traj_setup_EAM.xyz  \n" \
    f"run {step}  \n" \
    f"write_data {location}/data.halfhalf_EAM  \n" 
    with open(f'{location}/in.setup_pinning_auto_EAM' , "w") as f:
        f.write(param_string)
    return param_string


def write_in_pinning(location: str, pressure : float, temperature : float, step = 10000, hkl=[6,6,0],
IP_a=20.00, IP_kappa=4.0, barostat=1.00, thermostat=0.25):
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.halfhalf_EAM  \n" \
    "mass          1 48.0  \n" \
    "pair_style eam/fs  \n" \
    "pair_coeff * * Ti1.eam.fs Ti  \n" \
    "neighbor        0.3 bin \n" \
    "timestep        0.002 \n" \
    "run_style	verlet \n" \
    f"velocity all create  {temperature}  1 mom yes rot yes \n" \
    f"fix ensemble all npt temp {temperature} {temperature} {thermostat} z {pressure} {pressure} {barostat} \n" \
    "fix 100 all momentum 100 linear 1 1 1 \n" \
    f"fix bias all rhok {hkl[0]} {hkl[1]} {hkl[2]}   {IP_kappa}   {IP_a} \n" \
    "thermo_style custom step temp pzz pe lz f_bias f_bias[1] f_bias[2] f_bias[3] vol etotal \n" \
    "thermo 50 \n" \
    f"dump dumpXYZ all xyz 500 {location}/traj_pinning_EAM.xyz \n" \
    f"run {step} \n" 
    with open(f'{location}/in.pinning_auto_EAM' , "w") as f:
        f.write(param_string)
    return param_string


def write_in_setup_liquid(location: str, pressure : float, temperature : float, step = 10000, hkl=[6,6,0],
    barostat=1.00, thermostat=0.25):
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.crystal_end_meaned  \n" \
    "mass          1 48.0  \n" \
    "pair_style eam/fs  \n" \
    "pair_coeff * * Ti1.eam.fs Ti  \n" \
    "neighbor        0.3 bin \n" \
    "timestep        0.002 \n" \
    "run_style	verlet \n" \
    f"fix ensemble all nvt temp 4000.0 {temperature} {thermostat} \n" \
    f"fix bias all rhok {hkl[0]} {hkl[1]} {hkl[2]}  0.0 0.0 \n" \
    "thermo 50 \n" \
    "thermo_style custom step temp press density f_bias[3] vol etotal \n" \
    "#dump dumpXYZ all xyz 500 traj_liquid_melting_EAM.xyz \n" \
    f"run {step} \n" \
    f"write_data {location}/data.liquid_EAM \n"
    with open(f'{location}/in.setup_liquid_auto_EAM' , "w") as f:
        f.write(param_string)
    return param_string

def write_in_liquid(location: str, pressure : float, temperature : float, step = 10000, hkl=[6,6,0],
    barostat=1.00, thermostat=0.25):
    param_string="units           metal \n" \
    "dimension	3  \n" \
    "boundary        p p p  \n" \
    "atom_style  atomic  \n" \
    f"read_data     {location}/data.liquid_EAM  \n" \
    "mass          1 48.0  \n" \
    "pair_style eam/fs  \n" \
    "pair_coeff * * Ti1.eam.fs Ti  \n" \
    "neighbor        0.3 bin \n" \
    "timestep        0.002 \n" \
    "run_style	verlet \n" \
    f"velocity all create  {temperature} 1 mom yes rot yes \n" \
    f"fix ensemble all npt temp {temperature} {temperature} {thermostat} z {pressure} {pressure} {barostat} \n" \
    "fix 100 all momentum 100 linear 1 1 1 \n" \
    f"fix bias all rhok {hkl[0]} {hkl[1]} {hkl[2]}  0.0 0.0 \n" \
    "thermo_style custom step temp pzz pe lz f_bias f_bias[1] f_bias[2] f_bias[3] vol etotal enthalpy \n" \
    "thermo 50 \n" \
    f"dump dumpXYZ all xyz 500 {location}/traj_liquid_EAM.xyz \n" \
    f"run {step} \n" 
    with open(f'{location}/in.liquid_auto_EAM' , "w") as f:
        f.write(param_string)
    return param_string

class ReadLammps():
    #Read lammps output: assumes Step is outputted.
    def __init__(self, filein):
        self.get_dict(filein)
    def get_dict(self,filein):
        with open(filein) as f:
            lines=f.readlines()
        for i, val in enumerate(lines):
            if 'Step' in val:
                break
        start_loc=i+1        
        for i, val in enumerate(lines):
            if 'Loop time' in val:
                break
        end_loc=i

        labels = lines[start_loc-1].split()[0:]
        mydict = {}
        for i in labels:
            mydict[i]=[]
        for i, val in enumerate(lines[start_loc:end_loc]):
                for j, val2 in enumerate(labels):
                    mydict[val2].append(val.split()[j])
        self.dict = mydict
    def average_Q(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['f_bias[3]'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['f_bias[3]'], dtype=float).mean()

    def std_Q(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['f_bias[3]'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['f_bias[3]'], dtype=float).std()


    def average_vol(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Volume'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Volume'], dtype=float).mean()

    def average_density(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Density'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Density'], dtype=float).mean()

    def average_energy(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['TotEng'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['TotEng'], dtype=float).mean()

    def average_potential_pinned(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            total_pot=np.array(self.dict['PotEng'], dtype=float)[-thinned:]
            bias=np.array(self.dict['f_bias'], dtype=float)[-thinned:]
            return (total_pot - bias).mean()
        return (total_pot - bias).mean()

    def std_potential_pinned(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            total_pot=np.array(self.dict['PotEng'], dtype=float)[-thinned:]
            bias=np.array(self.dict['f_bias'], dtype=float)[-thinned:] #not sure if this is contained in PE or not...
            return (total_pot).std()
        return (total_pot - bias).std()

    def average_bias_pot(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['f_bias'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['f_bias'], dtype=float).mean()

    def std_bias_pot(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['f_bias'], dtype=float)[-thinned:].std()
        return np.array(self.dict['f_bias'], dtype=float).std()

    def average_temp(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Temp'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Temp'], dtype=float).mean()

    def average_pressure(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Press'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Press'], dtype=float).mean()

    def average_enthalpy(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Enthalpy'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Enthalpy'], dtype=float).mean()

    def average_lx(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Lx'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Lx'], dtype=float).mean()

    def average_ly(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Ly'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Ly'], dtype=float).mean()

    def average_lz(self, thinned=None):
        if thinned!=None and type(thinned) == int:
            return np.array(self.dict['Lz'], dtype=float)[-thinned:].mean()
        return np.array(self.dict['Lz'], dtype=float).mean()



class IPparams():
    def __init__(self, location, pressure, temperature, step, hkl, thermostat, barostat, z_half, Nz, a, kappa, thinned, N, atoms):
        self.location=location
        self.pressure=pressure
        self.barostat=barostat
        self.temperature=temperature
        self.thermostat=thermostat
        self.step=step
        self.hkl=hkl
        self.z_half=z_half
        self.Nz = Nz
        self.a=a
        self.kappa=kappa
        self.thinned=thinned
        self.N = N
        self.atoms = atoms
        return

    def write_all(self):
        self.write_crystal()
        self.write_liquid()
        self.write_coex()
        return
    
    def write_crystal(self):
        write_in_crystal(location=self.location, pressure=self.pressure, temperature=self.temperature, step=self.step,
        hkl=self.hkl, thermostat=self.thermostat, barostat=self.barostat)
        return

    def write_liquid(self):
        write_in_setup_liquid(location=self.location, pressure=self.pressure, temperature=self.temperature, step=self.step,
        hkl=self.hkl, thermostat=self.thermostat, barostat=self.barostat)
        write_in_liquid(location=self.location, pressure=self.pressure, temperature=self.temperature, step=self.step, hkl=self.hkl,
        thermostat=self.thermostat, barostat=self.barostat)
        return
    
    def write_coex(self):
        write_in_setup_pinning(location=self.location, pressure=self.pressure, temperature=self.temperature, step=self.step,
        hkl=self.hkl, z_half=self.z_half, thermostat=self.thermostat, barostat=self.barostat)
        write_in_pinning(location=self.location, pressure=self.pressure, temperature=self.temperature, step=self.step, hkl=self.hkl,
        IP_a=self.a, IP_kappa=self.kappa, thermostat=self.thermostat, barostat=self.barostat)
        return



    def update_IP(self, kappa_prefactor=1):
        crystal = ReadLammps(f'{self.location}/crystal_auto_EAM.out')
        liquid = ReadLammps(f'{self.location}/liquid_auto_EAM.out')
        q_crystal = crystal.average_Q(thinned=self.thinned)
        q_liquid = liquid.average_Q(thinned=self.thinned)
        self.a = q_liquid+ (q_crystal-q_liquid)/2
        self.kappa = kappa_prefactor*kB*self.temperature*self.Nz**2/(q_crystal-q_liquid)**2
        return



    
