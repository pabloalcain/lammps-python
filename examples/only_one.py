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

temp = np.linspace(2.0,0.5,31)
param = {}
param['T'] = 2.0
param['lambda'] = 20
param['N'] = 5488
param['x'] = 0.5
param['d'] = 0.04
param['potential'] = 'medium'
comp = ('rdf', 'ssf', 'mste', 'mink', 'thermo')

system = MDSys(silent=False)
system.setup(param, comp)
system.minimize()
#system.read_dump("dump.lammpstrj")
system.lmp.command("run 0")
for j in range(10):
  system.run(10)
  system.dump()
  system.results()
  print '\b\b\b\b\b{0:.>3}%'.format((j+1)*2),
  stdout.flush()
system.flush()
print '\b\b\b\b\bDone!\r'
