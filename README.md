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
Go to displit root directory and type :
```
$ cmake .
$ make
$ ./main X
```
where `X` is either 0 (open existing file and do statistics) or 1 (recompute 
simulation before doing statistics).

### References
[1] T. Liou, A.H. Mueller, S. Munier. 
Fluctuations of the multiplicity of produced particles in proton-nucleus collisions.
2016

[2] G.P. Salam.
OEDIPUS: Onium evolution, dipole interaction and perturbative unitarization simulation.
Comput.Phys.Commun. 105 (1997) 62-76.