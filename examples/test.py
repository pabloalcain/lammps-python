#!/usr/bin/python
"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
import numpy as np
from itertools import product

temp=np.linspace(1.6,0.1,6)
#temp=[]
T = 4.2
l = 4
N = 6
xarr = np.linspace(0.5,0.4,11)
darr = np.linspace(0.01,0.12,12)
V = "medium"
for x, d in product(xarr,darr):
        print "x={0}, d ={1}".format(x, d)
	system = MDSys(T, l, N, x, d, V, gpu=False)
	system.build_script(fname="test.inp")#, dump="esponja.lammpstrj")
	system.lmp.command("compute mste all mste/atom 5.4")
	system.setup(lind=False)

	for i in temp:
	    system.set_T(i)
	    for j in range(5):
		system.run(1)
		system.results(r_mink = 0.2, r_cell = 0.1)
		if ( j % 5 == 0): system.dump()
	    system.flush()
