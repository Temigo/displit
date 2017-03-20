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
*Warning* : filenames are used to get/set parameters. Do not change it.


#### Generation of events and histograms
###### generate [nb_events] [delta] [max_y] [cutoff_type]
Generate `nb_events` events with UV cutoff `delta`, maximal rapidity `max_y` and
an IR cutoff type `cutoff_type`. Each event is stored as a TTree in a ROOT file.
You can stop it with Ctrl-C and events generated until this point will be saved.

###### generate-mpi [parameters file]
Example of parameters file :
```
1000000 0.01 3.00 gaussian
1000000 0.01 3.00 lorentzian1
```
The pattern for each line is `nb_events delta max_y cutoff_type` where `delta` 
is the UV cutoff and `cutoff_type` is the IR cutoff.

###### fluctuations [filename]
###### fluctuations-mpi [file]
Example of files list (1 filename per line) :
```
mpi_tree_100000events_cutoff0.010000_ymax2.000000_gaussian.root
mpi_tree_100000events_cutoff0.010000_ymax2.000000_lorentzian1.root
```

###### generate-fluctuations-mpi [parameters file]
Example of parameters file :
```
100 10000 0.01 3.00 2.0 gaussian
100 10000 0.01 3.00 2.0 lorentzian1
```
The pattern for each line is `nb_tasks nb_events_per_task delta max_y R cutoff_type`.

#### Plots and graphics
###### draw-fluctuations [histogram file] [max_y]
Draw the distribution of probability of having n dipoles (fluctuations) after an evolution until
rapidity `max_y`. 

###### compare [histofile]
Draw the fluctuations for different histograms on the same canvas.

###### compare-c [histofile] [mode]
Draw the value of c for different values of parameter `mode` which can be `r`, 
`delta` or `ymax`. The histograms filenames are in `histofile`. All parameters
other than `mode` should be fixed.

###### compare-R [histofile]
Test the fluctuations dependance in x_{01}/R.

#### Miscellaneous
###### list-cutoffs
List available IR cutoffs.

###### draw-cutoffs
Draw some IR cutoffs (not all of them).

###### check [tree_file]
Check the real number of events recorded in [tree_file] against the theoretical number (encoded in
the filename).

###### fit-bare-r [delta] [max_y]
Draw the distribution of dipoles size for given parameters : UV cutoff `delta` and
maximal rapidity `max_y`, comparing with and without an IR cutoff (rigid).

###### ancestors [nb_events] [rho] [max_y]
Not fully tested yet.

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