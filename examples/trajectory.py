from QMDAnalysis.trajectories import TrajectorySH
from QMDAnalysis.molecules import Molecule
import matplotlib.pyplot as plt
import numpy as np

"""
Some examples of how the Trajectory class can be used to analyse a single surface hopping trajectory,
as well as perform scattering calculations and transformations.
"""

BOHR = 0.52917721092  # Angstroms


def internal_coordinates_example():
    """ example of how a trajectory can be read in and the internal coordinates calculated. """

    fpath = '../data/Trajectories/CS2/diss/TRAJ_00020/output.xyz'

    traj = TrajectorySH.init_from_xyz(fpath)
    print(f'Trajecectory from {fpath} has:')
    print(f'time steps: {traj.nts}, max time: {max(traj.time)}, dt: {traj.dt}')
    angles = traj.get_angles([1, 0, 2])
    r1 = traj.get_bond_lengths([1, 0])
    r1 *= BOHR
    r2 = traj.get_bond_lengths([2, 0])
    r2 *= BOHR # will be in a.u. by default

    fig, ax = plt.subplots(nrows=3, sharex=True)
    ax[0].plot(traj.time, angles)
    ax[0].set_ylabel('SCS angle (deg.)')
    ax[1].plot(traj.time, r1)
    ax[1].set_ylabel('CS distance (Ang.)')
    ax[2].plot(traj.time, r2)
    ax[2].set_ylabel('CS distance (Ang.)')
    ax[2].set_xlabel('t (fs)')
    fig.show()


def trajectory_rmsd_example():
    """ calculate the rmsd over time wrt a ref structure (here eq. CS2) """

    fpath = '../data/Trajectories/CS2/diss/TRAJ_00020/output.xyz'
    ref_path = '../data/Molecules/cs2.xyz'

    traj = TrajectorySH.init_from_xyz(fpath)
    ref = Molecule.init_from_xyz(ref_path)
    rmsd = traj.kabsch_rmsd(ref)

    fig, ax = plt.subplots()
    ax.plot(traj.time, rmsd)
    ax.set_ylabel('RMSD')
    ax.set_xlabel('t (fs)')
    fig.show()


def trajectory_electron_scattering():
    from QMDAnalysis.scattering import IAM_trajectory_scattering, IAM_form_factors

    fpath = '../data/Trajectories/CS2/diss/TRAJ_00010/output.xyz'
    traj = TrajectorySH.init_from_xyz(fpath)

    qmin, qmax = 0, 24
    Nq = 481
    qvec = np.linspace(qmin, qmax, Nq)

    mol = traj.geometries[0]
    FF, fq = IAM_form_factors(mol, qvec, ELEC=True) # XRS by default, if ELEC=True - electron scattering
    Itrj = IAM_trajectory_scattering(traj, qvec, fq, FF, ELEC=True)

    ref_signal = Itrj[0, :]

    # in x-ray case also divide through by ref_signal
    Idiff = (100 * 0.03) * (Itrj - ref_signal)

    tgrid, qgrid = np.meshgrid(traj.time, qvec)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(qgrid, tgrid, Idiff.T, cmap='RdBu')
    ax.set_xlabel('q (A^-1)')
    ax.set_xlim([0, 12])
    ax.set_ylabel('t (fs)')
    fig.show()


def trajectory_xray_scattering():
    from QMDAnalysis.scattering import IAM_trajectory_scattering, IAM_form_factors

    fpath = '../data/Trajectories/CS2/diss/TRAJ_00010/output.xyz'
    traj = TrajectorySH.init_from_xyz(fpath)

    qmin, qmax = 0, 24
    Nq = 481
    qvec = np.linspace(qmin, qmax, Nq)

    mol = traj.geometries[0]
    FF, fq = IAM_form_factors(mol, qvec)
    Itrj = IAM_trajectory_scattering(traj, qvec, fq, FF)

    ref_signal = Itrj[0, :]

    # in x-ray case also divide through by ref_signal
    Idiff = (100 * 0.03) * ((Itrj - ref_signal) / ref_signal)

    tgrid, qgrid = np.meshgrid(traj.time, qvec)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(qgrid, tgrid, Idiff.T, cmap='RdBu')
    ax.set_xlabel('q (A^-1)')
    ax.set_xlim([0, 12])
    ax.set_ylabel('t (fs)')
    fig.show()


if __name__ == "__main__":
    trajectory_electron_scattering()
