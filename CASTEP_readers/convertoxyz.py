import ase, ase.io
import sys
from warnings import simplefilter
simplefilter('ignore')

#example call: python3 converttoxyz.py *.castep
for input_raw in sys.argv[1:]:
    input_string=input_raw.split('.')[0]
    output_string = f'{input_string}.xyz'
    atoms=ase.io.read(input_raw)
    atoms.arrays.pop('initial_magmoms', None)
    atoms.arrays.pop('initial_charges', None)
    ase.io.write(output_string, atoms)
