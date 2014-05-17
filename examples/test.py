#!/usr/bin/python
"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
import numpy as np

temp=np.linspace(1.6,0.1,6)
#temp=[]
T = 4.2
l = 4
N = 6
x = 0.5
d = 0.2
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=False)
system.build_script(fname="test.inp")#, dump="esponja.lammpstrj")
system.lmp.command("compute mste all mste/atom 5.4")
system.setup(lind=False)

for i in temp:
    system.set_T(i)
    for j in range(50):
        system.run(100)
        system.results(r_mink = 0.2, r_cell = 0.1)
        system.log()
        system.dump()
    system.flush()
