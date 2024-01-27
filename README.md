# QMDAnalysis

A python module for analysis of mixed quantum-classical trajectory simulations, calculation of scattering observables, frequencey calculations, coordinate transforms, and the generation of potential energy surfaces.

## Installation

* git clone this repo
* run `pip install -e .` in the root directory of the repo for development
* otherwise run `pip install .` in the root directory

## Modules

### Trajectories:
* Initialise trajectories from xyz files
* Calculate internal coordinates
* Calculate RMSD over time wrt a reference structure

### Scattering:
* Calculate rotationally averaged scattering signals
* Non-resonant x-ray scattering
* Electron scattering

### Frequencies:
* Calculate vibrational frequencies and normal modes from a Hessian file
* Compute the Hessian of a structure through a wrapper to pyscf (dft, mp2, ccsd, scf)
* Displace a reference structure along selected normal modes

### Coordinate Transforms:
* Transform from cartesian to normal mode coordinates and back
* Transform from cartesian to internal coordinates

### Potential Generation:
* Generate grids that can be used in an electronic structure code
* Normal mode coordinates
* Cartesian coordinates
* Internal coordinates

### DVR:
* TBD


## Generating documentation

install pdoc3 if not installed:
```
pip install pdoc3
```

to generate a local server to view documentation of the modules run the following command in
the root directory of the repo
```
pdoc3 --http localhost:8080 QMDAnalysis 
```

go to:
```
http://localhost:8080
```
