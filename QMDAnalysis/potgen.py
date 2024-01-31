import numpy as np
from QMDAnalysis.molecules import Molecule, Frequency
from QMDAnalysis.constants import BOHR

class Potential:

    def __init__(self, eq_geom : Frequency):
        self.equilibrium = eq_geom
        self.cartesian_coordinates = None
        self.qcoordinates = None

    def generate_normal_mode(self, transformation_matrix, qmins, qmaxs, mode_indicies, ngrid):
        """
        Generate cartesian coordinates for electronic structure calculations from normal mode coordinates.
        The resulting grid is a direct product grid that is linear in normal mode coordinates.
        The total number of grid points is ngrid**n where n = len(mode_indicies).

        Parameters
        ----------
        transformation_matrix : numpy.ndarray
            transformation matrix from normal mode to cartesian coordinates, columns are eigenmodes.
        qmins: list
            minimum value of the normal mode coordinates to be varied - ordered according to mode_indicies
        qmaxs: list
            maximum value of the normal mode coordinates to be varied - ordered accoridng to mode_indicies
        mode_indicies: list
            which normal modes to vary. [0, 1, ... , n] are the first, second and nth normal mode according to
            energy ordering.
        ngrid: int
             number of points along each dimension in normal mode coordinates.

        """
        nvibs = self.equilibrium.nvibs
        ndof_total = self.equilibrium.natoms * 3
        ndof_rot_trans = ndof_total - nvibs
        ndof_active = len(qmins)
        if ndof_active != len(mode_indicies):
            raise Exception('Number of indices of variable vibrational modes must match the number of axis limits.')
        ngrid_total = ngrid**ndof_active

        grids = []
        for vib in range(ndof_active):
            qgrid = np.linspace(qmins[vib], qmaxs[vib], ngrid)
            grids.append(qgrid)
        meshed_grids = np.meshgrid(*grids, indexing='ij')
        grid_points = np.column_stack([axis.flatten() for axis in meshed_grids])
        qvib_coords = np.zeros((ngrid_total, nvibs))
        qvib_coords[:, mode_indicies] = grid_points
        qcoordinates = np.concatenate([np.zeros((ngrid_total, ndof_rot_trans)), qvib_coords], axis=1)

        cartesian_coords = np.zeros((self.equilibrium.natoms, 3, ngrid_total))
        for i in range(ngrid_total):
            cartesian_coords[:, :, i] = self._norm2cart(qcoordinates[i, :], transformation_matrix)

        self.qcoordinates = qcoordinates
        self.cartesian_coordinates = cartesian_coords

    def _norm2cart(self, qcoordinates, transformation_matrix):
        mw_dcoords = transformation_matrix @ qcoordinates
        umw_dcoords = mw_dcoords.reshape((self.equilibrium.natoms, 3)) * (1 / np.sqrt(self.equilibrium.Zs[:, np.newaxis]))
        cart_coords = umw_dcoords + self.equilibrium.coordinates
        return cart_coords

    def write_to_xyz(self, fname):
        natoms, _, npoints = self.cartesian_coordinates.shape
        with open(fname, 'a+') as f:
            for i in range(npoints):
                coords = self.cartesian_coordinates[:, :, i] * BOHR
                f.write(f'{natoms}\n')
                f.write(f'grid point {i+1}\n')
                array_string = '\n'.join(
                    [f'{self.equilibrium.atom_labels[ind]} ' + ' '.join(map(str, row)) for ind, row in enumerate(coords)]) + '\n'
                f.write(array_string)

