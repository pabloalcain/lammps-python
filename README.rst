LAMMPS Python wrapper
=====================

This project started [as can be seen from previous commits] as a
series of functions to execute LAMMPS in a python environment.
Initially, the idea was to be able to produce a thermalization scheme,
non-existent in LAMMPS. Eventually it grew, and now tries to focus in
being a LAMMPS wrapper. Because of the infrastructure and the kind of
work we do, it is meant to work on single core implementations with
GPU acceleration. This is obviously *against* LAMMPS idea of being
massively parallel, but take it as a proof of concept of the
interface. This is, still, a library interface meant to be used by our
group, with a bit of work to make it extendable. If this approach is
eventually desired as another mode of LAMMPS, a plugin, whatever, the
next obvious (and pretty hard) step will be to work with the typical
domain decomposition in LAMMPS. Probably, apart from the obvious and
already mentioned domain decomposition issue, there are other problems
due to the narrowness of the initial project.

Installation
------------

First build the analysis library:

.. code-block:: bash

   cd library/
   make

Then run, in the root directory

.. code-block:: bash

   python setup.py install

(root privileges needed)
