## DISPLIT Project
__Goal__ : Investigate hadrons internal structure through the multiplicity 
distribution of gluons in high-energy scatterings of hadrons.

Using Monte-Carlo simulations in C++/ROOT

### Install and configure ROOT framework
ROOT is a C++ framework for large scale data analysis developed at CERN.
See https://root.cern.ch/.

The CMakeLists.txt file assumes you have defined an environment variable called
`$ROOTSYS` pointing to your ROOT installation directory.

For instance you could add these lines in your `~/.bashrc` file:
```
export ROOTSYS=/usr/local
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
```

### Compile and run
You will also need MPI libraries to compile displit. For example on Ubuntu run :
```
$ sudo apt-get install libopenmpi-dev
```

Go then to displit root directory and type :
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./main [options]
```

To run displit, either use `./main [options]` or if using MPI 
(only with options `generate-mpi` and `fluctuations-mpi`) 
use `mpiexec -np $NB_TASKS ./main [options]`.

### Options
* generate [nb_events] [rho] [max_y] [cutoff_type]
* generate-mpi [parameters file]
Example of parameters file :
```
1000000 0.01 3.00 gaussian
1000000 0.01 3.00 lorentzian1
```
The pattern for each line is `nb_events rho max_y cutoff_type`.

* fluctuations [filename]
* fluctuations-mpi [file]
Example of files list :
```
mpi_tree_100000events_cutoff0.010000_ymax2.000000_gaussian.root
mpi_tree_100000events_cutoff0.010000_ymax2.000000_lorentzian1.root
```
Warning : the filename is used to get the parameters. Do not change it.

* draw-fluctuations [histogram file] [max_y]
* check [file]
* fit-bare-r [rho] [max_y]
* ancestors [nb_events] [rho] [max_y]
* draw-cutoffs

### Running displit with SLURM + MPI
The file `generic.slurm` is provided as an example to launch displit on a cluster
using SLURM and OpenMP.

### References
[1] T. Liou, A.H. Mueller, S. Munier. 
Fluctuations of the multiplicity of produced particles in proton-nucleus collisions.
2016

[2] G.P. Salam.
OEDIPUS: Onium evolution, dipole interaction and perturbative unitarization simulation.
Comput.Phys.Commun. 105 (1997) 62-76.