from QMDAnalysis.potgen import Potential
from QMDAnalysis.molecules import Frequency


def vary_bending_mode():
    hess_file = '../data/Freq/H2O_B3LPY_hessian.txt'
    coord_file = '../data/Molecules/h2o_b3lyp_opt.xyz'
    mol = Frequency.init_from_xyz(coord_file, units='angstroem')
    freqs, nmodes, transformation_mat = mol.calculate_normal_modes(hessian_file=hess_file)
    pot = Potential(mol)
    qmins = [-1.5]
    qmaxs = [1.5]
    variable_mode = [0]
    ngrid = 21
    pot.generate_normal_mode(transformation_mat, qmins, qmaxs, variable_mode, ngrid)
    pot.write_to_xyz('../data/Potential/h2o_nm0.xyz')


def vary_symm_mode():
    hess_file = '../data/Freq/H2O_B3LPY_hessian.txt'
    coord_file = '../data/Molecules/h2o_b3lyp_opt.xyz'
    mol = Frequency.init_from_xyz(coord_file, units='angstroem')
    freqs, nmodes, transformation_mat = mol.calculate_normal_modes(hessian_file=hess_file)
    pot = Potential(mol)
    qmins = [-2]
    qmaxs = [1]
    variable_mode = [1]
    ngrid = 21
    pot.generate_normal_mode(transformation_mat, qmins, qmaxs, variable_mode, ngrid)
    pot.write_to_xyz('../data/Potential/h2o_nm1.xyz')


def vary_asymm_mode():
    hess_file = '../data/Freq/H2O_B3LPY_hessian.txt'
    coord_file = '../data/Molecules/h2o_b3lyp_opt.xyz'
    mol = Frequency.init_from_xyz(coord_file, units='angstroem')
    freqs, nmodes, transformation_mat = mol.calculate_normal_modes(hessian_file=hess_file)
    pot = Potential(mol)
    qmins = [-1.5]
    qmaxs = [1.5]
    variable_mode = [2]
    ngrid = 21
    pot.generate_normal_mode(transformation_mat, qmins, qmaxs, variable_mode, ngrid)
    pot.write_to_xyz('../data/Potential/h2o_nm2.xyz')


def generate_whole_pes():
    hess_file = '../data/Freq/H2O_B3LPY_hessian.txt'
    coord_file = '../data/Molecules/h2o_b3lyp_opt.xyz'
    mol = Frequency.init_from_xyz(coord_file, units='angstroem')
    freqs, nmodes, transformation_mat = mol.calculate_normal_modes(hessian_file=hess_file)
    pot = Potential(mol)
    qmins = [-1.5, -2, -1.5]
    qmaxs = [1.5, 1, 1.5]
    variable_mode = [0, 1, 2]
    ngrid = 21
    pot.generate_normal_mode(transformation_mat, qmins, qmaxs, variable_mode, ngrid)
    pot.write_to_xyz('../data/Potential/h2o_pes.xyz')


if __name__ == "__main__":
    generate_whole_pes()