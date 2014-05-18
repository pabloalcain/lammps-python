#!/usr/bin/python
"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
import numpy as np
import sys
from sys import stdout


Ti = 1.6
Tf = 0.1
nT = 51
temp=np.linspace(Ti, Tf, nT)
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
    system.set_T(i, tdamp = 100.0)
    print "T: {0:>6}. Equilibrating...".format(i),
    stdout.flush()    
    system.equilibrate(wind = 20, nfreq = 1000)
    print "\bDone! Evolution.....0%",
    stdout.flush()
    for j in range(50):
        system.run(1000)
        system.results()
        system.dump()
        print '\b\b\b\b\b{0:.>3}%'.format((j+1)*2),
        system.flush()
    print '\b\b\b\b\bDone!'
