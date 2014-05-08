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

Install
-------

Requires numpy and lammps

To install, just add src/ to the environment variable PYTHONPATH

Analysis
========

In the folder analysis/, there are two routines separately that have the
RDF and Minkowski algorithms used. Lindemann is already to be made, and
MSTE has been added to LAMMPS code

Goal
----

This project just uses the LAMMPS library and implements some analysis
tools for post process. These should be later used with ctypes to bind 
it to Python INPUT process.

The routines that are going to be implemented here are:

- RDF: Radial distribution function, g(r)

- Minkowski: The 4 3-D Minkowski functionals

- Lindemann: Calculate Lindemann coefficient (not done)

The LAMMPS library is dynamically loaded, and in the simple.cpp file
we see how they are accessed.

How to use
----------

Change in the Makefile the include and library directories, and type

$ make

then run ./simple.e

to build as a shared library, just run

$ make library

