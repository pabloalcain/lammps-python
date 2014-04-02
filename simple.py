#!/usr/bin/python
#*-* coding: utf-8 *-*
import sys
sys.path.append("funciones")
import time
from os import chdir
from funciones import *
from funct_input import *
from numpy import linspace, logspace, log10
#import lammps
from lammps import lammps


argv = sys.argv
if len(argv) != 2:
  print "Syntax: simple.py in.lammps"
  sys.exit()

infile = sys.argv[1]
[inp,N,Ti,Tf,NTemp,sty,Nsteps,Ndump,fpath]=create_input(infile)

lmp = lammps("",["-c","off","-screen","data/lammps.out"])

# run infile one line at a time
lines = inp.splitlines()
for line in lines:
  lmp.command(line)

if (sty==1): #lineal
  temparr=linspace(Ti,Tf,NTemp)
elif (sty==2):#log
  temparr=logspace(log10(Ti),log10(Tf),NTemp)
elif (sty==3):#gauss
  temparr=gauss_samp(Ti,Tf,0.4,0.75,NTemp-1)
#  temparr=gauss_samp(Ti,Tf,0.3,0.5,NTemp-1)
else:
  raise ValueError("Nunca deberíamos entrar acá")

for i in temparr:
    set_temp(lmp,i,fpath)
    Therm(lmp,i,fpath,N)
    Run(lmp,Nsteps,Ndump,i,fpath)
