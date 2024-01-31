#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --time=01:00:00

module purge
module load GCC/11.2.0 OpenMPI/4.1.1
module load ORCA/5.0.3

cp $SLURM_SUBMIT_DIR/*.inp /tmp
cp $SLURM_SUBMIT_DIR/*.xyz /tmp

cd /tmp && $EBROOTORCA/bin/orca orca_calculation.inp > orca_calculation.out

cp /tmp/*.gbw $SLURM_SUBMIT_DIR
cp /tmp/*.xyz $SLURM_SUBMIT_DIR
cp /tmp/*.out $SLURM_SUBMIT_DIR
