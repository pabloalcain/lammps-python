#!/usr/bin/python
"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
import numpy as np

temp=np.linspace(1.6,0.1,2)
#temp=[]
T = 1.6
l = 20
N = 5488
x = 0.5
d = 0.05
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=True)
system.build_script(fname="test.inp")#, dump="esponja.lammpstrj")
system.lmp.command("compute mste all mste/atom 5.4")
system.setup(lind=False)

for i in temp:
    system.set_T(i)
    system.equilibrate(wind = 10, nfreq = 100)
    for j in range(5):
        system.run(100)
        system.results()
        system.dump()
    system.flush()
