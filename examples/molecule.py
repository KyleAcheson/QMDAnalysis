from QMDAnalysis.molecules import Molecule
import numpy as np

""" Some basic examples of how the Molecule class can be used. """

# Build a molecule from scratch
labels = ['C', 'O', 'O']
coordinates = np.array([[0, 0, 0], [0, 0, 1.162], [0, 0, -1.162]])
m = Molecule(labels, coordinates) # if no units are provided - Angstroem is assumed

# Read a molecule from xyz file
fpath = '../data/Molecules/methane.xyz'
m = Molecule.init_from_xyz(fpath)

# calculate centre of mass and subtract from coordinates
com = m.centre_of_mass()
m.coordinates -= com

# Read from a molden file
fpath = '../data/Freq/cs2.molden'
m = Molecule.init_from_molden(fpath)

bl = m.bond_length([0, 1])
print(f'The C-S bond length (Bohr) is: {bl}')
ba = m.angle([1, 0, 2])
print(f'The S-C-S bond angle is: {ba}')

distance_matrix = m.distance_matrix()
print('The distance matrix (Bohr):')
print(distance_matrix)

fpath2 = '../data/Molecules/cs2_n2.xyz'
m1 = Molecule.init_from_xyz(fpath2)

print('The structure of molecule m: ')
for atom in m:
    print(atom)

print('The structure of molecule m1: ')
for atom in m1:
    print(atom)

# calculate the rmsd between two slightly different structures
rmsd = Molecule.kabsh_rmsd(m, m1)
print(f'The rmsd is between m and m1 is: {rmsd}')