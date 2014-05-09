from ctypes import *
import numpy as np
from lammps import lammps

analysis = cdll.LoadLibrary("libanalysis.so")

def minkowski(lmp, rad, rcell):
    tmp = (c_double * 4)()
    ir = c_double(rad)
    ic = c_double(rcell)
    analysis.minkowski(lmp, ir, ic, tmp)
    return np.fromiter(tmp, dtype = np.float, count = 4)


def rdf(lmp, nbins, rmax):
    npair = 4
    tmp = (c_double * (nbins * (1 + npair * 2)))()
    inb = c_int(nbins)
    ir = c_double(rmax)
    analysis.rdf(lmp, inb, ir, tmp)
    for i in range(100):
      print i, 
      for j in range(9):
	print tmp[i*9+j],
      print ""


    print "haoa@"
    return np.frombuffer(tmp, dtype = np.float, count = nbins * npair)
