import ase, ase.io
import sys
from warnings import simplefilter
from ase import Atoms
simplefilter('ignore')


label=str(sys.argv[1])
force_name='QM_forces'
energy_name='QM_energy'
virial_name='QM_virial'
#example call: python3 bulkcasteptoxyz.py label_name dir/*.castep
atoms=[]
for i, input_raw in enumerate(sys.argv[2:]):
    input_string=input_raw.split('.')[0]
    atoms.append(ase.io.read(input_raw))

for i, a in enumerate(atoms):
    a.arrays[f'{force_name}'] = a.get_forces()
    a.info[f'{energy_name}'] = a.get_potential_energy(force_consistent=True)
    a.info[f'{virial_name}'] = -a.get_stress(voigt=False)*a.get_volume()
    a.info[f'config_type'] = f'{label}'

Properties=f'species:S:1:pos:R:3:{force_name}:R:3'
ase.io.write(f'{label}.extxyz', atoms, columns=['symbols', 'positions', f'{force_name}'], format='extxyz')