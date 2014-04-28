"""
This file continues a previous dump file
"""
from md_classes import MDSys
from lammps import lammps
from numpy import linspace

T = 0.6
l = 20
N = 5000
x = 0.5
d = 0.05
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=True)
system.build_script(fname="lammps.inp", dump="minim.lammpstrj")
system.setup()
system.set_T(0.6, therm = True) 
system.run(50000)
#  system.results()
# RUN Nsteps steps
# POSTPROCESS
