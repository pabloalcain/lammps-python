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
from sys import stdout

param = {}
param['lambda'] = 20
param['N'] = 60
param['x'] = 0.4
param['potential'] = 'medium'
temp = np.linspace(2.0,0.5,31)
darr = np.linspace(0.01,0.08,8)
comp = ('rdf', 'ssf', 'mink', 'mste', 'thermo', 'fit')
for d in darr:
    param['d'] = d
    param['T'] = 2.0
    print "x={0} d={1}".format(0.4, d)
    system = MDSys(silent=False, mste = "mste" in comp)
    system.setup(param, comp)
    system.minimize()
    for T in temp:
        param['T'] = T
        system.setup(param, comp)
        system.dump()
        print "T: {0:>6}. Equilibrating...".format(T),
        stdout.flush()    
        system.equilibrate(wind = 20, nfreq = 1000)
        print "\bDone! Evolution.....0%",
        stdout.flush()
        for j in range(1):
            system.run(1000)
            system.results()
            system.dump()
            print '\b\b\b\b\b{0:.>3}%'.format((j+1)*2),
            stdout.flush()
        system.flush()
        print '\b\b\b\b\bDone!\r'
