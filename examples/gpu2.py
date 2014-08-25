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
from itertools import product
from time import time
from datetime import timedelta

T = 0.5
l = 20
N = 5488
x = 0.5
d = 0.05
V = 'medium'
system = MDSys(T, l, N, x, d, V, gpu=True)
system.build_script(fname="test.inp", dump="gpu2.lammpstrj")
system.lmp.command("compute mste all mste/atom 5.4")
system.setup(lind=False, path='./med-5r')

system.set_T(T, tdamp = 100.0)
system.run(10)
system.results()
system.flush()
