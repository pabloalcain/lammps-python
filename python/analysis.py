"""
This is simply a wrapper for the library libanalysis.so
"""

from ctypes import *
import numpy as np
#from lammps import lammps

analysis = cdll.LoadLibrary("libanalysis.so")

def minkowski(lmp, rad, rcell):
    """
    Wrapper for minkowski. Need to add a warning with the
    size of the lattice for digitalization!
    """

    tmp = (c_double * 4)()
    analysis.minkowski.argtypes = [c_void_p, c_double, c_double, c_void_p]
    analysis.minkowski(lmp, rad, rcell, tmp)
    return np.fromiter(tmp, dtype = np.float, count = 4)


def rdf(lmp, nbins, rmax, npairs):
    """
    Wrapper for the c++ libanalysis. Its main disadvantage is that we
    should know beforehand the number of pairs to consider.  Since
    this is meant to Neutron Stars, we document here behavior when
    npair is set to 4:

    If we label the pairs with index k, so that
    k = 1 => all v all
    k = 2 =>  1  v  1
    k = 3 =>  1  v  2
    k = 4 =>  2  v  2

    Returns a 2D array, with this format.

    1st column: Distance
    1 + 2 * k: g(r) between the pair labeled by k
    2 + 2 * k: int(g(r)) between the pair labeled by k
    """

    ncol = npairs * 2 + 1
    tmp = ( c_double * (nbins * ncol) ) ()
    analysis.rdf.argtypes = [c_void_p, c_int, c_double, c_void_p]
    analysis.rdf(lmp, nbins, rmax, tmp)
    r = np.frombuffer(tmp, dtype = np.float, count = nbins * ncol)
    return np.reshape(r, (nbins, ncol))
