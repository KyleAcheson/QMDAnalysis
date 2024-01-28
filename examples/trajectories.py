from QMDAnalysis.trajectories import Ensemble
from QMDAnalysis.molecules import Molecule
import numpy as np
import os
import matplotlib.pyplot as plt

def ensemble_internal_coordinates():

    """ Load an ensemble from file and calculate the bond angles of each trajectory
        as well as the average bond angle over time. """

    parent_dir = '../data/Trajectories/CS2/diss'
    traj_paths = [subdir.path + '/output.xyz' for subdir in os.scandir(parent_dir) if subdir.is_dir()]

    ensemble = Ensemble.load_ensemble(traj_paths, traj_type='sh')
    print(f'Loaded {ensemble.ntrajs} trajectories into ensemble')

    bond_angles = np.full((ensemble.ntrajs, ensemble.nts_max), np.nan)
    for i, traj in enumerate(ensemble):
        print(f'trajectory {i+1} ran for: {traj.nts} time steps')
        bond_angle = traj.get_angles([1, 0, 2])
        bond_angles[i, :traj.nts] = bond_angle
        fig, ax = plt.subplots()
        ax.plot(traj.time, bond_angles[i, :traj.nts])
        ax.set_xlabel('t (fs)')
        ax.set_ylabel('SCS angle (deg.)')
        ax.set_title(f'trajectory {i+1}')
        fig.show()

    avg_angles = np.nansum(bond_angles, axis=0)
    avg_angles /= ensemble.tcount
    fig, ax = plt.subplots()
    ax.plot(np.linspace(0, ensemble.nts_max*0.5, ensemble.nts_max), avg_angles)
    ax.set_xlabel('t (fs)')
    ax.set_ylabel('SCS angle (deg.)')
    ax.set_title('average trajectory')
    fig.show()


def ensemble_rmsd_example():
    parent_dir = '../data/Trajectories/CS2/bound'
    ref_file = '../data/Molecules/cs2.xyz'
    traj_paths = [subdir.path + '/output.xyz' for subdir in os.scandir(parent_dir) if subdir.is_dir()]
    ensemble = Ensemble.load_ensemble(traj_paths, traj_type='sh')
    print(f'Loaded {ensemble.ntrajs} trajectories into ensemble')

    ref_mol = Molecule.init_from_xyz(ref_file)
    rmsd = ensemble.kabsch_rmsd(ref_mol)

    fig, ax = plt.subplots()
    for i in range(ensemble.ntrajs):
        ax.plot(rmsd[i, :], label=f'traj {i+1}')
    ax.set_ylabel('RMSD')
    ax.set_xlabel('time step')
    ax.set_title('RMSD wrt eq. structure')
    ax.legend(frameon=False)
    fig.show()

    avg_rmsd = np.sum(rmsd, axis=0) / ensemble.tcount
    fig, ax = plt.subplots()
    ax.plot(avg_rmsd)
    ax.set_ylabel('RMSD')
    ax.set_xlabel('time step')
    ax.set_title('Average ensemble RMSD')
    fig.show()

    breakpoint()



if __name__ == "__main__":
    ensemble_rmsd_example()