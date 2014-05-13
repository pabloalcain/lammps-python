"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from md_classes import MDSys
from lammps import lammps
from numpy import linspace


temp=linspace(1.6,0.1,6)
temp=[4.2]
T = 4.2
l = 20
N = 5488
x = 0.5
d = 0.05
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=False)
system.build_script(fname="test.inp", dump="esponja.lammpstrj")
system.lmp.command("compute mste all mste/atom 5.4")
system.setup()
for i in temp:
    system.set_T(i, therm = False) 
    system.run(0)
    mste = system.mste()
    m = system.minkowski(5.5, 2.1)
    r = system.rdf(200, 23.94)
    th = system.thermo()
    s = system.structure(r)
    #  system.results()
# RUN Nsteps steps
#Step Temp KinEng E_pair TotEng Press
#0   0.50174298   0.75247734    15.156092    15.908569   0.88371423    
