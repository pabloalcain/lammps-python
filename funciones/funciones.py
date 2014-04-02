#*-* coding: utf-8 *-*
from numpy import *
import random
import os

def gauss(x,s,m):
  y=zeros((size(x),1))
  if size(x)==1:
    if x==0:
      return 0
    else:
      return sqrt(1/(2*pi*s**2))*exp(-(x-m)**2/(2*s**2))
  for i in range(size(x)):
    if x[i]==0:
      y[i]=0
    else:
      y[i]=sqrt(1.0/(2*pi*s**2))*exp(-(x[i]-m)**2/(2*s**2))
  return y


def cgauss(x,s,m):
  y=zeros((size(x),1))  
  for i in range(size(x)):
    if i == 0:
      y[i]=gauss(x[i],s,m)
    else:
      y[i]=y[i-1]+gauss(x[i],s,m)
  return y/y[size(x)-1]

def gauss_samp(ti,tf,s,m,Nsteps):
  x=linspace(ti,tf,1000)
  a=cgauss(x,s,m)
  l=[]
  for i in range(Nsteps):
    for j in range(size(a)):
      if a[j]>float(i)/Nsteps:
        l.append(x[j])
        break
  l.append(x[-1])
  return l 

def Run(lmp,Nsteps,Ndump,temp,fpath):
    print "#Corriendo evolución de %i pasos a temperatura %g" % (Nsteps, temp)
    lmp.command("reset_timestep 0")
    lmp.command("log %sT%g/evol.log" % (fpath,temp))
    lmp.command("dump\t\t1 all custom %i %sT%g/evol.lammpstrj type id x y z vx vy vz" % (Ndump, fpath,temp))
    lmp.command("dump_modify\t1 sort id")
    lmp.command("run %i" % Nsteps)
    lmp.command("undump 1")
    print "#Completado!"
    
def Therm(lmp,temp,fpath,Npart):
    st_break=10*temp**2/Npart
    print "#Termalizando a temperatura %g, objetivo %g" % (temp,st_break)
    lmp.command("reset_timestep 0")
    lmp.command("log %sT%g/therm.log" % (fpath,temp))
    lmp.command("dump\t\t1 all custom %i %sT%g/therm.lammpstrj type id x y z vx vy vz" % (100000, fpath, temp))
    lmp.command("dump_modify\t1 sort id")
    N=100
    alpha=0.001
    while 1:
        std=0.0
	sxy=0.0
	sx=0.0
	sy=0.0
	sxx=0.0
        for i in range(N):
            lmp.command("run 1000")
            t=lmp.extract_compute("thermo_temp",0,0)
	    e=lmp.extract_compute("thermo_pe",0,0)
	    print i, e/Npart, t, temp
	    sxy+=e/Npart*i
	    sx+=i
	    sy+=e/Npart
	    sxx+=i**2
            std+=(t-temp)**2
        std=std/N
	k=(N*sxy-sx*sy)/(N*sxx-sx**2)
	print "#std=%g,k=%g"% (std,k)
        if (std<st_break and k>-alpha*abs(sx)/N**2): break
    lmp.command("undump 1")
    print "#Completado!"

def set_temp(lmp,temp,fpath):
    seed=random.randint(1,10000)
    os.makedirs(fpath+"T%g"%temp)
    #lmp.command("fix 2 all langevin %g %g 100.0 %d" % (temp, temp, seed))
    lmp.command("fix 2 all nvt temp %g %g 10.0" % (temp, temp))
