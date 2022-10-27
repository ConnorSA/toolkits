import ase
import numpy as np
class AtomsCustom(ase.Atoms):
    """"ase.Atoms object for implementation of custom methods"""
    def split_atom(self, atom_index: int, species: list, seperation = None, direction = None):
        """Split selected atom into two, returns new ase.Atoms object.
        
        atom_index : atom to split.
        species    : ['X', 'Y'] where X and Y are the chemical species of the new atoms
        seperation : distance between atoms X and Y, default is found as bond seperation between nearest neighbour.
        direction  : [x,y,z] vector to move atoms by +/- seperation/2 centred on atom that was split.
        """
        import numpy as np
        if seperation == None:
            self.seperation = np.min(self.get_distances(0, range(1, len(self.get_positions())), mic=True))
        else: self.seperation=seperation

        if direction !=None:
            direction = np.array(direction)
            if len(direction) == 3: self.direction = direction/np.linalg.norm(direction) 
            else: print('direction = [x,y,z]')            
        else:
            direction = np.random.rand(3) - 0.5
            direction = direction/np.linalg.norm(direction)
            self.direction = direction
        centroid = self.get_positions()[atom_index]
        x0 = centroid + self.direction*self.seperation/2
        x1 = centroid - self.direction*self.seperation/2
        split_atoms=self.copy()
        split_atoms.pop(atom_index)
        split_atoms = split_atoms + ase.Atoms(species,positions=[x0,x1],cell=self.get_cell())
        return 
        
def one_interstitial_fcc(atoms : ase.atoms, sc : list):
    """Return a list containing fcc crystal with 1 interstitial
    atoms : ASE atoms object that is the FCC orthorhombic (conventional unit cell).
    sc    : Size of the supercell, sc = [nx,ny,nz], where ni is the number of repeating unit cells.
    """
    interstitials=[]
    sc = np.array(sc)
    atoms_fcc_sc=atoms*sc
    isolated_atom=atoms.copy()
    for i in range(len(isolated_atom)-1):
        isolated_atom.pop()
    isolated_atom.set_cell(atoms_fcc_sc.cell)
    octahedral_site=np.array([0.5,0.5,0.5])/sc
    isolated_atom.set_scaled_positions([octahedral_site])
    interstitials.append(atoms_fcc_sc+isolated_atom)
    tetrahedral_site=np.array([0.25,0.25,0.25])/sc
    isolated_atom.set_scaled_positions([tetrahedral_site])
    interstitials.append(atoms_fcc_sc+isolated_atom)
    return interstitials
