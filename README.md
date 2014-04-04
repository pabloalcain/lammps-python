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