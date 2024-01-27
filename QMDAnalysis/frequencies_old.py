import numpy as np
from itertools import permutations
from molecules import Molecule
from constants import *


class Frequencey(Molecule):

    def __init__(self, labels, coordinates, units='angstroem'):
        super().__init__(labels, coordinates, units)
        if self._is_linear():
            self.nvibs = 3 * self.natoms - 5
        else:
            self.nvibs = 3 * self.natoms - 6

    def _is_linear(self):
        if self.natoms > 3:
            return False
        elif self.natoms == 2:
            return True
        else:
            elements = [i for i in range(self.natoms)]
            atom_permutations = permutations(elements)
            for indicies in atom_permutations:
                angle = self.angle(indicies)
                if abs(180 - angle) <= 1E-5:
                    return True
            return False

    def calculate_normal_modes(self, **kwargs):
        if kwargs.get('hessian_file'):
            hessian_file = kwargs['hessian_file']
            hessian = np.genfromtxt(hessian_file)
            freqs, normal_modes, transformation_matrix = self._normal_mode_calculator(hessian)
        else:
            hessian = self._pyscf_hessian(**kwargs)
            freqs, normal_modes, transformation_matrix = self._normal_mode_calculator(hessian)
        return freqs, normal_modes, transformation_matrix

    def _normal_mode_calculator(self, hessian):

        fconstants_au, modes = self.diag_hessian(hessian)
        freqs = np.sqrt(np.abs(fconstants_au))  # in a.u.
        freqs_wavenums = freqs * AU2Hz / LIGHT_SPEED_SI * 1e-2
        normal_modes = np.einsum('z,zri->izr', self.Zs ** -.5, modes.reshape(self.natoms, 3, -1))
        transformation_matrix = modes.T
        return freqs_wavenums[-self.nvibs:], normal_modes, transformation_matrix

    def _construct_mass_matrix(self):
        masses = np.array(self.Zs)
        mass_vec = np.repeat(masses, 3)
        mass_mat = np.sqrt(np.outer(mass_vec, mass_vec))
        return mass_mat

    def diag_hessian(self, hessian):
        sq_mass_matrix = self._construct_mass_matrix()
        weighted_hessian = hessian * (1 / sq_mass_matrix)
        fconstants, modes = np.linalg.eigh(weighted_hessian)
        return fconstants, modes

    def _pyscf_hessian(self, **kwargs):
        from pyscf import gto

        method = kwargs.get('method')
        if not method:
            raise Exception('A method (dft, scf, mp2, ccsd) must be provided.')

        pyscf_geom = self.coordinates.tolist()
        for i in range(self.natoms):
            pyscf_geom[i].insert(0, self.atom_labels[i])

        mol = gto.M(atom=pyscf_geom, basis=kwargs['basis'], verbose=0)

        if method == 'dft':
            hessian = self._call_dft(mol, kwargs['xc'])
        elif method == 'scf':
            hessian = self._call_scf(mol)
        elif method == 'ccsd':
            hessian = self._call_ccsd(mol)
        elif method == 'mp2':
            hessian = self._call_mp2(mol)
        else:
            raise Exception(f'{method} is not currently implemented.')

        hessian = hessian.transpose(0, 2, 1, 3).reshape(self.natoms*3, self.natoms*3)
        return hessian

    @staticmethod
    def _call_scf(mol):
        from pyscf import scf, hessian
        mf = scf.RHF(mol)
        mf.kernel()
        return mf.Hessian().kernel()

    @staticmethod
    def _call_dft(mol, xc):
        from pyscf import dft, hessian
        mf = dft.RKS(mol)
        mf.xc = xc
        mf.kernel()
        hessian = mf.Hessian().kernel()
        return hessian

    @staticmethod
    def _call_mp2(mol):
        from pyscf import scf, mp, hessian
        rhf = scf.RHF(mol)
        rhf.kernel()
        mp2 = mp.MP2(rhf)
        return mp2.kernel()

    @staticmethod
    def _call_ccsd(mol):
        from pyscf import scf, cc, hessian
        rhf = scf.RHF(mol)
        rhf.kernel()
        ccsd = cc.CCSD(rhf)
        return ccsd.kernel()

    


if __name__ == "__main__":
    hessian_file = '../data/H2O_B3LPY_hessian.txt'
    fname = '../data/eq_geom.xyz'
    ref_geom = Frequencey.init_from_xyz(fname)
    freqs, normal_modes, tmat = ref_geom.calculate_normal_modes(hessian_file=hessian_file)
    #freqs, normal_modes, tmat = ref_geom.calculate_normal_modes(method='dft', xc='B3LYP', basis='def2-svp')
    breakpoint()