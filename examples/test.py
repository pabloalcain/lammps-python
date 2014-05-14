"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from md_classes import MDSys
from lammps import lammps

from numpy import linspace
from numpy import concatenate, savetxt
from numpy import array


temp=linspace(1.6,0.1,6)
temp=[4.2]
T = 4.2
l = 5
N = 10
x = 0.5
d = 0.2
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=False)
system.build_script(fname="test.inp")#, dump="esponja.lammpstrj")
system.lmp.command("compute mste all mste/atom 5.4")
system.setup()

fout=open("prueba",'w')

for i in temp:
    system.set_T(i, therm = False) 
    headers="#vol, sur, br, eul, temp, ke, epair, etot, press"
    print>>fout,headers
    for j in range(50):
        system.run(1000)
        t_c = system.mste()
        [vol, sur, br, eul] = system.minkowski(0.2, 0.1)
        t_r = system.rdf(200, 23.94)
        [temp, ke, epair, etot, press] = system.thermo()
        t_s = system.structure(t_r)

        print>>fout,vol, sur, br, eul, temp, ke, epair, etot, press
        system.dump()
        
    #  system.results()
