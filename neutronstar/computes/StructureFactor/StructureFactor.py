"""
Compute class
"""

from neutronstar.computes import Compute
import ctypes as ct
import numpy as np
import os
_DIRNAME = os.path.dirname(__file__)
libssf = ct.CDLL(os.path.join(_DIRNAME, 'libssf.so'))
ssf_c = libssf.ssf

class StructureFactor(Compute):
  """
  Structure Factor computation from definition.
  """
  def __init__(self, k, pairs, repetitions=2, lebedev=194):
    """
    Calculate the structure factor from angular definition,
    integrating with Lebedev quadrature

    Parameters
    ----------

    k : iterable (preferrably numpy array)
        Wavenumbers to calculate

    pairs: iterable
        List of pairs to consider. Each element of the list is a tuple
        of size 2, with the first element being the types considered
        and the second one the other types. For example, if we want to
        calculate the RDF on a system with 2 types of particles of all
        types vs all types and type 1 vs type 2, the pairs should be:

        pairs = [((1, 2), (1, 2)), ((1,), (2,))]

        A keyword for "all types" is type 0


    repetitions : int, optional
        Number of repetitions of the principal cell to
        consider. Default value is 2

    lebedev : int, optional
        Number of points in the sphere in which to calculate the
        Lebedev quadrature. Default value is 194

    Notes
    -----

    Points in the sphere cannot take any value, but only previously
    tabulated ones.

    Pairs have to be of same types, we cannot calculate by definition
    (yet, at least) the structure factor of 1 - 2 types.

    So far we cannot put more than one key inside each tuple, since it
    would probably either some serious refactoring or messy code. 'All
    particles' can still be chosen by keyword value 0
    """
    self.pairs = pairs
    self.rep = repetitions
    self.k = k
    self.lebedev = lebedev
    super(StructureFactor, self).__init__()
    self.header = []
    for p in pairs:
      h = []
      for t in p:
        h.append(','.join(map(str, t)))
      self.header.append('-'.join(map(str, h)))

  def compute(self, system):
    """
    Calculate structure factor.

    Parameters
    ----------

    system : System
        System on which we calculate the Structure Factor

    Returns
    -------

    ssf : dict
        An array with the information of the compute. The first column
        is 'r' and the rest is the structure factor calculated for the
        pair list.
    """
    natoms = system['N']
    size = system.size
    npairs = len(self.pairs)
    ncol = npairs + 1
    pair_ar = np.zeros(2*npairs)
    i = 0
    for p in self.pairs:
      pair_ar[i] = p[0][0]
      pair_ar[i+1] = p[1][0]
      i += 2
    ncol = len(self.pairs) + 1
    npoints = len(self.k)
    tmp = (ct.c_double * (npoints * ncol))()
    x_p = system.x.ctypes.data_as(ct.c_void_p)
    t_p = system.t.ctypes.data_as(ct.c_void_p)
    k_p = self.k.ctypes.data_as(ct.c_void_p)
    pair_p = pair_ar.ctypes.data_as(ct.c_void_p)
    ssf_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int,
                          ct.c_double, ct.c_int, ct.c_int, ct.c_int,
                          ct.c_void_p, ct.c_int, ct.c_void_p,
                          ct.c_void_p]
    ssf_c(x_p, t_p, natoms, size, npoints, self.lebedev, self.rep,
              k_p, tmp)
    ssf = np.frombuffer(tmp, dtype=np.double, count=npoints * ncol)
    return ssf.reshape((npoints, ncol))
