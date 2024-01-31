Contains coordinates for generating the ground state
potential energy surface of H2O at the B3LYP/def2-SVP level.

The script generate_inputs.sh is an example of how
to generate all the input electronic structure
calculations from the list of xyz geometries in h2o_pes.xyz.

This shell script takes the xyz file, an Orca template file,
and a SLURM submission script as input.

NOTE: This is just an example, and only works for a specific
SLURM setup using Orca. Modify it to work on your HPC config
and with the electronic structure code of choice.
