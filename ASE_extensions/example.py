from ASE_ext import AtomsCustom
import ase, ase.build
Ti= ase.build.bulk('Ti')
Ti_sc = Ti*[4,1,1]
atoms_custom = AtomsCustom(Ti_sc)
atoms_split=atoms_custom.split_atom(0, ['Ti','Al'])
print(atoms_split)
