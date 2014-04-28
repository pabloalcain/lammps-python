"""
Example run in neutron stars with gpu
"""
from md_classes import MDSys
from lammps import lammps
from numpy import linspace

temp=linspace(1.6,0.1,60)
T = 1.6
l = 20
N = 5000
x = 0.5
d = 0.05
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=True)
system.build_script(fname="lammps.inp")
system.setup()
for i in temp:
    system.set_T(i, therm = True) 
    system.run(50000)
#  system.results()
# RUN Nsteps steps
# POSTPROCESS
