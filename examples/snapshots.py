"""
Example run in neutron stars with gpu.
Create snapshots 
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
import numpy as np
from itertools import product
from sys import stdout, argv

param = {}
param['potential'] = 'medium'
param['lambda'] = 20
param['N'] = 5488
param['x'] = float(argv[1])
darr = np.linspace(0.01,0.08,8)
temp = np.linspace(2.0,0.5,31)
comp = ('rdf', 'ssf', 'mink', 'mste', 'thermo', 'fit')

for d in darr:
  param['d'] = d
  print "x={0} d={1}".format(x, d)
  system = MDSys(gpu=True, root="./snapshots")
  system.setup(param, comp)
  system.minimize()
  for T in temp:
    param['T'] = T
    system.setup(param, comp)
    print "\rT: {0:>6}. Equilibrating...".format(T),
    stdout.flush()
    system.equilibrate(wind = 20, nfreq = 1000)
    system.dump()
    print 'Done!'
