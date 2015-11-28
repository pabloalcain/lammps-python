Analysis
========

It's a tiny library that performs some typical analysis. It comes bundled
with an executable that shows how to use them to interact with LAMMPS

The routines that are going to be implemented here are:

- RDF: Radial distribution function, g(r)

- Minkowski: The 4 3-D Minkowski functionals

- MSTE: The clusterization algorithm with energy binding. (not done/moved to LAMMPS)

- Lindemann: Calculate Lindemann coefficient (not done)

How to use
----------

Change in the Makefile the include and library directories, and type

```make executable```

then run `./simple.e`

to build as a shared library, just run

```make library```

and to install the library,

```make install```
