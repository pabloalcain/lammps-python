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

It can be discussed how this postprocess should be done. The first idea
is to mimic the old behaviour: LAMMPS create the dumpfiles, and then it
is post-processed, with a different class/routine to get:

- RDF: Radial distribution function, g(r), and its FFT, S(k)

- Minkowski: The 4 3-D Minkowski functionals

- MSTE: The clusterization algorithm with energy binding.

- Lindemann: Calculate Lindemann coefficient

Eventually, we could put the bindings *inside* the run (design should be
discussed and thought) or, even better, as LAMMPS computes
