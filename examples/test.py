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
N = 6
x = 0.5
d = 0.0005
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=False)
system.build_script(fname="test.inp")#, dump="esponja.lammpstrj")
system.lmp.command("compute mste all mste/atom 5.4")
system.setup()
for i in temp:
    system.set_T(i, therm = False) 
    for j in range(50):
        system.run(1000)
        mste = system.mste()
        m = system.minkowski(5.5, 2.1)
        r = system.rdf(200, 23.94)
        th = system.thermo()
        s = system.structure(r)
        system.dump()
        
#        system.log()
    #  system.results()
# RUN Nsteps steps
#real    0m5.490s
#user    0m4.905s
#sys     0m0.077s
