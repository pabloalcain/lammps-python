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
    tmp = POINTER(c_double * nbins)
    inb = c_int(nbins)
    ir = c_double(rmax)
    analysis.rdf(lmp, inb, ir, tmp)
    return np.fromiter(tmp, dtype = np.int, count = nbins)
