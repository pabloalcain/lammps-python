"""
This benchmark is to check initalization time
and penalty of splitting multiple runs
"""
from neutronstar.MDSys import MDSys
from lammps import lammps
from numpy import linspace
import time

def bench(tot, nsteps):
    t1 = time.time()
    for i in xrange(tot/nsteps):
        system.run(nsteps)
    t2 = time.time()
    print>>f, '{tot}, {n}, {TPS}'.format(tot = tot,
					 n = nsteps,
					 TPS = float(tot)/(t2-t1))

f = open('info', 'w')
T = 1.6
l = 20
N = 5000
x = 0.5
d = 0.05
V = "medium"
system = MDSys(T, l, N, x, d, V, gpu=True)
system.build_script(fname="lammps.inp")
system.setup()

print>>f, '#TotalSteps, StepsPerRun, TPS(TimestepsPerSecond)'

tot = 50000
steps = [100, 500, 1000, 2000, 5000, 10000, 50000]
for i in steps: bench(tot, i)
