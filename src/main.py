"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from md_classes import MDSys
from lammps import lammps
from numpy import linspace

temp=linspace(1.6,0.1,60)
T = 1.7
l = 20
N = 5000
x = 0.5
d = 0.05
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=False)
system.build_script(fname="lammps.inp")
system.setup()
for i in temp:
    system.set_T(i, therm = True) 
    system.run(100)
#  system.results()
# RUN Nsteps steps
# POSTPROCESS
