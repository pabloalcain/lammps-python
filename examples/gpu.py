"""
Example run in neutron stars with gpu
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
from numpy import linspace

temp=linspace(1.6,0.1,51)
T = 1.6
l = 20
N = 5000
x = 0.5
d = 0.8
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=True)
system.build_script(fname="gpu.inp")
system.setup()
for i in temp:
    system.set_T(i, therm = True) 
    system.run(50000)
    system.results(rmax = 80)
    system.dump()
    system.flush()
# RUN Nsteps steps
# POSTPROCESS
