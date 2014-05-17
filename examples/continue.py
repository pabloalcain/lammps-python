"""
This file continues a previous dump file
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
from numpy import linspace

T = 0.6
l = 20
N = 5000
x = 0.5
d = 0.05
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=True)
system.build_script(fname="continue.inp", dump="continue.lammpstrj")
system.setup(lind = False, mste = False)
system.set_T(0.6, therm = True) 

for j in range(50):
    system.run(1000)
    system.results()
    system.log
    system.dump()
system.flush()

#  system.results()
# RUN Nsteps steps
# POSTPROCESS
