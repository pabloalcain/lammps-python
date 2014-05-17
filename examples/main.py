"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
from numpy import linspace


temp=linspace(1.6,0.1,6)
T = 4.2
l = 20
N = 5000
x = 0.5
d = 0.005
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=False)
system.build_script(fname="main.inp")
system.lmp.command("compute mste all mste/atom 5.4")
system.setup(lind = False)
for i in temp:
    system.set_T(i)
    system.run(100)
    system.results()
    system.log()
    system.dump()
    system.flush()
    #  system.results()
# RUN Nsteps steps
