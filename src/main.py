"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from md_classes import MDSys
from lammps import lammps

T = 0.001
l = 20
N = 4
x = 0.5
d = 0.04
V = "medium"
system = MDSys(T, l, N, x, d, V)
system.build_script(fname="lammps.inp")
system.setup()
system.run(50000)

# RUN Nsteps steps
# POSTPROCESS
