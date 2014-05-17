Nuetron Star Input
==================

Goal
----

The idea of this project is to create a good framework to call LAMMPS
from Python, within the Neutron Star typical simulations. The core idea
is that we have a 6 dimensional parameters space, which are:

- V: Type of potential (Medium, Stiff, Horowitz)

- l: lambda, screening length of the Coulomb-Debye interaction

- N: Number of particles

- x: fraction of protons

- d: number density of particles

- T: temperature

It will have *no* input file, at least in its core. Python flexibility
allows an easy extension to support this behaviour, but the idea is to
create python files that instantiates the classes that will be developed
in this code.

Postprocess
-----------

It can be discussed how this postprocess should be done. Nevertheless,
so far we have an analysis library that is dynamically loaded and wrapped
in analysis.py. We have:


- RDF: Radial distribution function, g(r), and its FFT, S(k)

- Minkowski: The 4 3-D Minkowski functionals

- MSTE: The clusterization algorithm with energy binding.

- Lindemann: Calculate Lindemann coefficient (not done yet)

Install
-------

First install the analysis library:

$ cd src/
$ make install

and make sure libanalysis.so is in your path

Then run, in the root directory

$ python setup.py install

(root privileges needed)

Warning
-------

In order to use compute mste/atom, you should have built
LAMMPS with the compute_mste_atom.cpp/h files. These files
are in the lammps folder
