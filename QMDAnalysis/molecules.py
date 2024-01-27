from typing import Tuple
import numpy as np
import numpy.typing as npt
from copy import deepcopy
from itertools import permutations

from QMDAnalysis.constants import *

__all__ = [
    'Molecule',
    'Frequency'
]


class MoleculeTypeError(TypeError):
    def __init__(self, msg='Requires a Molecule object', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

class AngleDefError(ValueError):
    def __init__(self, msg='Angles are defined between three atoms', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

class DihedralDefError(ValueError):
    def __init__(self, msg='Dihedrals are defined between two bond planes (4 atoms)', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

class XYZTypeError(TypeError):
    def __init__(self, msg='File type must be in .xyz format', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

class MoldenTypeError(TypeError):
    def __init__(self, msg='File type must be in .molden format', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)


class Molecule:
    """
    A static molecule object. Includes all information on a molecular structure at a given point
    in time along a trajectory. Currently does not include any values of momentum, but can be extended to include.
    The structure can be built in a script by the user using its constructor which requires a list of atomic labels
    for each atom and a set of molecular coordinates stored as a numpy.ndarray where the rows corrospond to the atom number.
    Alternatively (as is the norm), one can construct the object using the class constructors `init_from_xyz`
    and `init_from_molden` which reads in coordinates and atomic labels from an external xyz or molden file.

    One can calculate all, or a selection of the internal coordinates for the molecule using the internal coordinate methods.
    This class also includes methods for calculating the rmsd between two molecule structures via. the Kabsch algorithm.
    It is also possible to calculate observables for the static molecule, although currently only using the scattering
    module. With the addition of fitting modules, it will be possible to fit an observable from a single static
    structure to an experimental observable.

    Attributes
    ----------
    atom_labels : list
        a list of strings containing atomic symbols for each atom in the molecule - forced to be upper case.
    coordinates : numpy.ndarray
        molecule cartesian coordinates - rows corrosond to the atoms and columns xyz.
    natoms : int
        total number of atoms in molecule.
    nelecs : list
        list of the number of electrons for each atom entry of the coordinate matrix.
    Zs : list
        list of the atomic masses for each atom entry of the coordinate matrix .
    """

    def __init__(self, labels, coordinates, units='angstroem'):
        self._atom_labels = self._check_atom_labels(labels)
        self.natoms, self.nelecs, self.Zs = self._get_atoms_info()
        if units == 'angstroem':
            self._coordinates = coordinates * (1 / BOHR)
        else:
            self._coordinates = coordinates

    def __repr__(self):
        return f"Molecule({self._atom_labels}, {self.natoms}, {self.nelecs}, {self.Zs}, {self._coordinates.__repr__()})"

    def __iter__(self):
        return MoleculeIterator(self)


    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coords: npt.NDArray):
        """
        Coordinate property is allowed to be set incase the instantiated coordinates are transformed.
        This setter ensures transformed coordinates are of same dimensionality.
        """
        self._check_array_coords(self.natoms, coords, 'Molecular coordinates')
        self._coordinates = coords

    @property
    def atom_labels(self):
        return self._atom_labels

    @staticmethod
    def _check_atom_labels(labels: list):
        """
        Checks atomic labels are correct type and upper case
        """
        if len(labels) != len(list(filter(str, labels))):
            raise Exception("Atomic labels must be a list of strings for each element")
        labels = [label.upper() for label in labels]
        return labels

    def _get_atoms_info(self) -> tuple[int, npt.NDArray, npt.NDArray]:
        """
        A private method to get atomic information
        :return number atoms, number electrons per atom, mass of each atom:
        :rtype: int, list, list
        """
        nelec_list, Z_list = [], []
        natoms = len(self._atom_labels)
        for atom in self._atom_labels:
            Z_list.append(periodic_table[atom]['Z'])
            nelec_list.append(periodic_table[atom]['nelec'])
        return (natoms, np.array(nelec_list), np.array(Z_list))

    @staticmethod
    def _check_array_coords(natoms: int, array: npt.NDArray, property: str):
        """
        Private method to check that property is defined for all molecular coordinates
        and has the correct type and dimensionality.
        """
        if type(array) != np.ndarray:
            raise Exception("%s must be specified as a numpy array" % property)
        dims = np.shape(array)
        if len(dims) != 2 and dims[1] != 3:
            raise Exception("%s must be an array with dimensions (natom, 3)" % property)
        if dims[0] != natoms:
            raise Exception("%s must be defined for %d atoms." % (property, natoms))

    @classmethod
    def init_from_xyz(cls, fpath: str, units='angstroem'):
        return cls._init_from_xyz(fpath, units)


    @classmethod
    def _init_from_xyz(cls, fpath: str, units='angstroem'):
        """
        A class constructor that instantiates a Molecule object from an external xyz file specified in the path.

        Parameters
        ----------
        fpath : str
            absolute path to the xyz file containing the molecular structure
        """
        atom_names, coords = cls.read_xyz_mol(fpath)
        return cls(
            labels=atom_names,
            coordinates=coords,
            units=units
        )

    @staticmethod
    def read_xyz_mol(fname: str) -> tuple[list, npt.NDArray]:
        """
        A static method which reads coordinates from an xyz file.

        Parameters
        ----------
        fname : str
            path to file.

        Returns
        -------
        labels: list
            a list of atomic labels for each atom in the coordinate matrix.
        atom_coords: numpy.ndarray
            an array that corrosponds to cartesian coordinates, rows are the atoms.

        """
        skip_lines = 1
        ftype = fname.split('.')[-1]
        if ftype != 'xyz':
            raise XYZTypeError('Molecule.init_from_xyz() requires a file in .xyz format')
        atom_coords, labels = [], []
        with open(fname, 'r') as f:
            first_line = f.readline().strip()
            Nat = int(first_line)
            for idx, line in enumerate(f):
                if idx >= skip_lines and idx <= Nat+skip_lines:
                    atom_list = line.strip().split()
                    atom_coords.append([float(i) for i in atom_list[1:]])
                    labels.append(atom_list[0])
        atom_coords = np.array(atom_coords)
        return labels, atom_coords

    @classmethod
    def init_from_molden(cls, fpath: str, units='angstroem'):
        return cls._init_from_molden(fpath, units)

    @classmethod
    def _init_from_molden(cls, fpath: str, units='angstroem'):
        """
        A class constructor that instantiates a Molecule object from an external molden file.

        Parameters
        ----------
        fpath : str
            path to the molden file
        """
        atom_labels, coords = cls.read_molden(fpath)
        return cls(
            labels=atom_labels,
            coordinates=coords,
            units=units
        )

    @staticmethod
    def read_molden(fname: str) -> tuple[list, npt.NDArray]:
        """
        A static method which reads coordinates a molden file.

        Parameters
        ----------
        fname : str
            path to file.

        Returns
        -------
        labels: list
            a list of atomic labels for each atom in the coordinate matrix.
        atom_coords: numpy.ndarray
            an array that corrosponds to cartesian coordinates, rows are the atoms.
        """
        ftype = fname.split('.')[-1]
        if ftype != 'molden':
            raise MoldenTypeError('Molecule.read_molden() requires a file in .molden format')

        labels, atoms = [], []
        mfile = open(fname, 'r')
        Atoms = False
        for line in mfile:
            if '[' in line or '--' in line:
                Atoms = False
            if '[Atoms]' in line:
                Atoms = True
            elif Atoms:
                words = line.split()
                labels += words[0]
                atom_vec = words[3:6]
                atoms += [[eval(coords) for coords in atom_vec]]

        atoms = np.array(atoms)
        return labels, atoms

    def distance_matrix(self):
        """
        A method to calculate the distance matrix of a given structure.

        Returns
        -------
        dist_mat : np.NdArray
            upper triangular distance matrix stored in sparse format (see scipy docs).
            To convert to a numpy.ndarray use the method `dist_mat.toarray()`.
        """
        dist_mat = np.zeros((self.natoms, self.natoms), dtype=np.float64)
        for i in range(self.natoms):
            for j in range(i + 1, self.natoms):
                rvec = self.coordinates[i, :] - self.coordinates[j, :]
                d = np.linalg.norm(rvec)
                dist_mat[i, j] = d
                dist_mat[j, i] = d
        return dist_mat

    def bond_length(self, bond_connectivity: list) -> float:
        """
        A method to calculate the bond length between a pair of atoms in the self.

        Parameters
        ----------
        bond_connectivity : list
            a list that contains the indexes of the pairs of atoms in the coordinate matrix.
            `len(bond_connectivity) == 2`

        Returns
        -------
        bond_len : float
            value of the bond length in angstrom
        """
        rvec = self.coordinates[bond_connectivity[0], :] - self.coordinates[bond_connectivity[1], :]
        bond_len = np.linalg.norm(rvec)
        return bond_len

    def angle(self, angle_connectivity: list) -> float:
        """
        A method to calculate the angle between a pair of bond lengths R_ij and R_kj.

        Parameters
        ----------
        angle_connectivity : list
            list that contains the indexes of the atoms that make up the angle.
            `len(angle_connectivity) = 3`, where the second index is the central atom i.e. j

        Returns
        -------
        theta : float
            the bond angle in degrees
        """

        if len(angle_connectivity) != 3 or self.natoms < 3:
            raise AngleDefError
        i, j, k = angle_connectivity
        r_ij = self.coordinates[i, :] - self.coordinates[j, :]
        r_kj = self.coordinates[k, :] - self.coordinates[j, :]
        cosine_theta = np.dot(r_ij, r_kj)
        sin_theta = np.linalg.norm(np.cross(r_ij, r_kj))
        theta = np.arctan2(sin_theta, cosine_theta)
        theta = 180.0 * theta / np.pi
        return theta

    def dihedral(self, dihedral_connectivity: list) -> float:
        """
        A method to calculate the dihedral angle between two bond lengths that form a plane.

        Parameters
        ----------
        dihedral_connectivity : list
            connectivity of the atoms that form the dihedral angle (i, j, k, l).
            `len(dihedral_connectivity) = 4`

        Returns
        -------
        phi: float
            dihedral angle in degrees

        """
        if len(dihedral_connectivity) != 4 or self.natoms < 4:
            raise DihedralDefError
        i, j, k, l = dihedral_connectivity
        r_ji = self.coordinates[j, :] - self.coordinates[i, :]
        r_kj = self.coordinates[k, :] - self.coordinates[j, :]
        r_lk = self.coordinates[l, :] - self.coordinates[k, :]
        v1 = np.cross(r_ji, r_kj)
        v1 /= np.linalg.norm(v1)
        v2 = np.cross(r_lk, r_kj)
        v2 /= np.linalg.norm(v2)
        p1 = np.cross(v1, r_kj) / np.linalg.norm(r_kj)
        a = np.dot(v1, v2)
        b = np.dot(p1, v2)
        phi = np.arctan2(b, a)
        phi = -180.0 - 180.0 * phi / np.pi
        if phi < -180.0:
            phi += 360.0
        return phi

    def centre_of_mass(self):
        """
        A function to calculate the centre of mass in cartesian space. Take away the resulting vector
        from a molecule to shift the origin to the centre of mass.

        Returns
        -------
        centre_mass: numpy.ndarray
            centre of mass in cartesian coordinates
        """
        tot = np.zeros((1, 3))
        for i in range(self.natoms):
            tot = tot + self.Zs[i] * self.coordinates[i, :]
        centre_mass = tot / np.sum(self.Zs)
        return centre_mass

    @staticmethod
    def kabsh_rmsd(molecule, referance_structure, Hydrogens=True, Mirror=False):
        """
        A method to calculate the minimum RMSD between two geometries through the Kabsch algorithm.
        This works by calculating the centroid of each vector X (i.e. `sum(x)/ len(x)`) and aligning the two
        geometries. Then by calculating the covariance matrix of the two centred structures, this is used
        to calculate the rotation matrix that minimises the rmsd through a procedure based on single value decomposition.

        Parameters
        ----------
        referance_structure : Molecule
            the seconnd molecule to with which the RMSD is calculated wrt the instance of the molecule objecy

        Returns
        -------
        lrms: float
            the lowest possible RMSD between the structures (after rotation and translation)
        """

        if not isinstance(referance_structure, Molecule):
            raise MoleculeTypeError('Kabsch algorithm requires the reference to be another molecular structure')
        if molecule.natoms != referance_structure.natoms:
            raise MoleculeTypeError('The two molecules must have the same dimensions')

        # new_molecule = deepcopy(molecule)
        # new_referance_structure = deepcopy(referance_structure)

        if not Hydrogens:
            molecule, referance_structure = Molecule.__remove_hydrogens_pair(molecule, referance_structure)

        aligned_molecule, aligned_ref = molecule.centre_align(referance_structure)
        lrms = aligned_molecule._lrmsd(aligned_ref)

        if Mirror:
            molecule.reflect(axis=0)
            aligned_molecule, aligned_ref = molecule.centre_align(referance_structure)
            lrms_mirrored = aligned_molecule._lrmsd(aligned_ref)
            molecule.reflect(axis=0)
            if lrms_mirrored < lrms:
                # print("Mirrored isomers detected")
                lrms = lrms_mirrored

        return lrms

    def _lrmsd(self, reference_molecule):
        diff = self.coordinates - reference_molecule.coordinates
        lrms = float(np.sqrt((np.sum(diff ** 2)) / self.natoms))
        return lrms

    def centre_align(self, referance_structure):
        """
        Takes an instance of Molecule and given some reference Molecule translate the coordinates to align centriod with
        the origin and find the optimal rotation matrix that aligns the Molecule instance with respect to reference_structure.

        Parameters
        ----------
        referance_structure : Molecule
            reference structure to align the Molecule instance with

        Returns
        -------
        molecule: Molecule
            a new instance of the original Molecule which has been aligned at the centroid and rotated to overlap the reference geometry
        ref: Molecule
            a new instance of the input reference_structure that has been aligned and rotated
        """
        molecule, ref = deepcopy(self), deepcopy(referance_structure)
        nc = np.shape(molecule.coordinates)[1]
        p0 = np.sum(molecule.coordinates, axis=0) / self.natoms
        q0 = np.sum(ref.coordinates, axis=0) / ref.natoms
        molecule.coordinates = molecule.coordinates - p0
        ref.coordinates = ref.coordinates - q0  # translate coords to align centroid w origin
        cov = np.transpose(molecule.coordinates) @ ref.coordinates  # calculate covariance matrix
        v, s, wh = np.linalg.svd(cov)  # do single value decomp. on covariance matrix
        w = wh.T
        w = np.squeeze(w)
        v = np.squeeze(v)
        eye = np.eye(nc)  # init identity matrix
        if np.linalg.det(w @ np.transpose(v)) < 0:
            eye[nc - 1, nc - 1] = -1
        u = w @ eye @ np.transpose(v)  # rotation matrix that minimises the rmsd
        for i in range(molecule.natoms):
            molecule.coordinates[i, :] = u @ molecule.coordinates[i, :]
        return molecule, ref

    def reflect(self, axis):
        if axis == 0:
            self.coordinates[:, 0] = -self.coordinates[:, 0]
        elif axis == 1:
            self.coordinates[:, 1] = -self.coordinates[:, 1]
        elif axis == 2:
            self.coordinates[:, 2] = -self.coordinates[:, 2]
        else:
            raise IndexError("Only three axes present, choose one of x, y or z")

    def remove_hydrogens(self):
        """
        Takes an instance of Molecule and returns a new Molecule instance which is a copy of the
        original instance but with the hydrogens removed.

        Returns
        -------
        mol_noh: Molecule
            copy of Molecule instance with the hydrogens removed
        """
        new_mol, new_labels = [], []
        for i in range(self.natoms):
            if self.atom_labels[i] != 'H':
                new_mol.append(self.coordinates[i, :])
                new_labels.append(self.atom_labels[i])
            else:
                pass
        mol_noh = Molecule(new_labels, np.array(new_mol))
        return mol_noh

    @staticmethod
    def __remove_hydrogens_pair(mol1, mol2):
        """Removes hydrogens for two molecule pairs - only used internally in kabsh_rmsd"""
        mol_a, mol_b, labels_a, labels_b = [], [], [], []
        for i in range(mol1.natoms):
            if mol1.atom_labels[i] != 'H':
                mol_a.append(mol1.coordinates[i, :])
                labels_a.append(mol1.atom_labels[i])
            if mol2.atom_labels[i] != 'H':
                mol_b.append(mol2.coordinates[i, :])
                labels_b.append(mol2.atom_labels[i])
            else:
                pass
        new_mol_a = Molecule(labels_a, np.array(mol_a))
        new_mol_b = Molecule(labels_b, np.array(mol_b))
        return new_mol_a, new_mol_b


class MoleculeIterator:

    def __init__(self, molecule):
        self._molecule = molecule
        self._index = 0

    def __next__(self):
        if self._index < len(self._molecule.coordinates):
            result = (self._molecule.atom_labels[self._index], self._molecule.coordinates[self._index, :])
            self._index += 1
            return result
        else:
            raise StopIteration
        
        

class Frequency(Molecule):

    def __init__(self, labels, coordinates, units='angstroem'):
        super().__init__(labels, coordinates, units)
        if self._is_linear():
            self.nvibs = 3 * self.natoms - 5
        else:
            self.nvibs = 3 * self.natoms - 6

    def freqs_from_molden(self, fname):
        molden_out = self._read_molden_nmodes(fname)
        freqs, normal_modes = molden_out[2], molden_out[3]
        return freqs, normal_modes

    @staticmethod
    def _read_molden_nmodes(fname: str) -> tuple[list, npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        A method to read coordinates, frequencies and normal modes from a molden files.
        Parameters
        ----------
        fname : str
            absolute path to molden file

        Returns
        -------
        labels : list
            a list of atom labels for each entry in the reference structures cartesian coordinate matrix
        atoms : numpy.ndarray
            array of coordinates for each atom in the reference structure
        freqs : numpy.ndarray
            the values of the frequencies of each mode in cm^-1
        vibs : numpy.ndarray
            the normal modes of each of the frequencies in contained in `freqs`
        """
        ftype = fname.split('.')[-1]
        if ftype != 'molden':
            raise MoldenTypeError('Vibration object must be instantiated from a .molden type file')

        labels, atoms, freqs, vibs = [], [], [], []

        mfile = open(fname, 'r')
        Atoms = False
        FREQ = False
        FRNORMCOORD = False

        actvib = -1
        for line in mfile:
            # what section are we in
            if '[' in line or '--' in line:
                Atoms = False
                FREQ = False
                FRNORMCOORD = False

            if '[Atoms]' in line:
                Atoms = True
            elif '[FREQ]' in line:
                FREQ = True
            elif '[FR-NORM-COORD]' in line:
                FRNORMCOORD = True
                # extract the information in that section
            elif Atoms:
                words = line.split()
                labels += words[0]
                atom_vec = words[3:6]
                atoms += [[eval(coords) for coords in atom_vec]]
            elif FREQ:
                freqs += [eval(line)]
            elif FRNORMCOORD:
                if 'vibration' in line or 'Vibration' in line:
                    vib_list = []
                    actvib += 1
                    if actvib > -1:
                        vibs += [vib_list]
                else:
                    vib_list += [[eval(coor) for coor in line.split()]]

        freqs = np.array(freqs)
        vibs = np.array(vibs)
        atoms = np.array(atoms)
        return (labels, atoms, freqs, vibs)

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
        """

        Calculates normal modes and their frequencies (in cm^-1).
        The calculation can be performed using a hessian read from file using the keyword argument *hessian_file*.
        Alternatively, the hessian can be calculated through a call to pyscf using the keyword arguments *method* and *basis*.
        If *method='dft'*, the keyword argument *xc* is also required.

        Parameters
        ----------
        kwargs

        hessian_file : str
            path to a file containing the hessian - normal modes calculated internally
        method : str
            level of theory (dft, mp2, ccsd, scf) at which to calculate the hessian at if not provided
        basis : str
            basis set used in hessian calculation
        xc : str
            exchange-correlation functional (for dft only)


        Returns
        -------

        """
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
        """
        Diagonalize the Hessian (will be mass weighted).

        Parameters
        ----------
        hessian : npt.NDArray
            unweighted hessian matrix

        Returns
        -------
        fconstants: npt.NDArray
            force constants in a.u.
        modes: npt.NDArray
            eigenvectors of the mass weighted Hessian (stored as columns)
        """
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
    
    @staticmethod
    def freq2time(freq_wavenum: float):
        """
        Convert a vibrational freq given in cm^-1 to a time period in fs

        Parameters
        ----------
        freq_wavenum : float
            vibrational freq in inverse cm

        Returns
        -------
        time_period: float
            period of the vibration in femtoseconds

        """
        lambd = 0.01 /freq_wavenum # wavelength in m
        time_period = 1 / (LIGHT_SPEED_SI /lambd ) # period of vibration in time (s)
        return time_period*1e15 # convert to fs



if __name__ == "__main__":
    fname = '../data/Molecules/h2o_b3lyp_opt.xyz'
    mol = Molecule.init_from_xyz(fname, 'angstroem')
    mol.angle([1, 0, 2])