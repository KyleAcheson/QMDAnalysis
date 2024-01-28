from QMDAnalysis.molecules import Frequency
import numpy as np


def example1():
    # build a molecule from scratch
    labels = ['C', 'O', 'O']
    coordinates = np.array([[0, 0, 0], [0, 0, 1.162], [0, 0, -1.162]])
    mol = Frequency(labels, coordinates, units='angstroem') # converts to bohr internally

    # calculate normal modes and frequencies through a call to pyscf - DOES NOT OPTIMIZE THE STRUCTURE!
    freqs, nmodes, transformation_mat = mol.calculate_normal_modes(method='dft', basis='def2-svp', xc='B3LYP')
    print('frequencies (cm^-1):')
    print(freqs)
    print('vibrational period (fs):')
    for freq in freqs:
        print(mol.freq2time(freq))


def example2():
    hess_file = '../data/Freq/H2O_B3LPY_hessian.txt'
    coord_file = '../data/Molecules/h2o_b3lyp_opt.xyz'
    mol = Frequency.init_from_xyz(coord_file, units='angstroem')
    freqs, nmodes, transformation_mat = mol.calculate_normal_modes(hessian_file=hess_file)
    print('frequencies (cm^-1):')
    print(freqs)
    print('vibrational period (fs):')
    for freq in freqs:
        print(mol.freq2time(freq))


def example3():
    # same as example2 but calculating the Hessian with pyscf instead of from file
    coord_file = '../data/Molecules/h2o_b3lyp_opt.xyz'
    mol = Frequency.init_from_xyz(coord_file, units='angstroem')
    freqs, nmodes, transformation_mat = mol.calculate_normal_modes(method='dft', basis='def2-svp', xc='B3LYP')
    print('frequencies (cm^-1):')
    print(freqs)
    print('vibrational period (fs):')
    for freq in freqs:
        print(mol.freq2time(freq))


if __name__ == "__main__":
    print('Building a CO2 calculation using pyscf')
    example1()
    print('H2O calculation using Hessian from file')
    example2()
    print('H2O calculation using Hessian calculated by pyscf')
    example3()
