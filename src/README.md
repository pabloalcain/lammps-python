Analysis
========

Goal
----

This project just uses the LAMMPS library and implements some analysis
tools for post process. These should be later used with ctypes to bind 
it to Python INPUT process.

The routines that are going to be implemented here are:

- RDF: Radial distribution function, g(r)

- Minkowski: The 4 3-D Minkowski functionals

- MSTE: The clusterization algorithm with energy binding. (not done)

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
