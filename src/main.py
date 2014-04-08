"""
This file is an example of how the MDSys class should be imported and
used. The main idea is that everything can be called like this, and
therefore looped from within Python
"""
from md_classes import MDSys
from lammps import lammps

T = 1.0
l = 20
N = 4000
x = 0.5
d = 0.04
V = "medium"
system = MDSys(T, l, N, x, d, V)
system.build_script(fname="lammps.inp")
lmp = lammps("",["-echo", "screen"])
system.setup(lmp)
lmp.command("run 1000")

# RUN Nsteps steps
# POSTPROCESS
log		data/{V}/l{l}/x{x}/N{N}/d{d}/T{T}/thermo.log
dump		1 all custom {ndump} &
		data/V{V}/l{l}/x{x}/N{N}/d{d}/T{T}/dump.lammpstrj &
		type id x y z vx vy vz
dump_modify	1 sort id
