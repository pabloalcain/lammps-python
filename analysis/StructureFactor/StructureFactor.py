"""
Calculation of the structure factor
"""

import ctypes as ct
import numpy as np
import os

_DIRNAME = os.path.dirname(__file__)
libssf = ct.CDLL(os.path.join(_DIRNAME, 'libssf.so'))
ssf_c = libssf.ssf

def structureFactor(x, t, size, pairs, k, rep=2, lebedev=194):
  """
  Calculate structure factor.

  Parameters
  ----------

  x : numpy array
      Positions of the particles in the system

  t : numpy array
      Types of the particles

  size : float
      Length of the box

  pairs: iterable
      List of pairs to consider. Each element of the list is a tuple
      of size 2, with the first element being the types considered
      and the second one the other types. For example, if we want to
      calculate the RDF on a system with 2 types of particles of all
      types vs all types and type 1 vs type 2, the pairs should be:

      pairs = [((1, 2), (1, 2)), ((1,), (2,))]

      A keyword for "all types" is type 0

  k : numpy array
      Wavenumbers to calculate

  repetitions : int, optional
      Number of repetitions of the principal cell to
      consider. Default value is 2

  lebedev : int, optional
      Number of points in the sphere in which to calculate the
      Lebedev quadrature. Default value is 194

  Returns
  -------

  ssf : dict
      An array with the information of the compute. The first column
      is 'r' and the rest is the structure factor calculated for the
      pair list.
  """
  natoms = np.shape(x)[0]
  npairs = len(pairs)
  ncol = npairs + 1
  pair_ar = np.zeros(2*npairs)
  i = 0
  for p in pairs:
    pair_ar[i] = p[0][0]
    pair_ar[i+1] = p[1][0]
    i += 2
  npoints = len(k)
  tmp = (ct.c_double * (npoints * ncol))()
  x_p = x.ctypes.data_as(ct.c_void_p)
  t_p = t.ctypes.data_as(ct.c_void_p)
  k_p = k.ctypes.data_as(ct.c_void_p)
  pair_p = pair_ar.ctypes.data_as(ct.c_void_p)
  ssf_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int,
                    ct.c_double, ct.c_int, ct.c_int, ct.c_int,
                    ct.c_void_p, ct.c_int, ct.c_void_p,
                    ct.c_void_p]
  ssf_c(x_p, t_p, natoms, size, npoints, lebedev, rep,
        pair_p, npairs, k_p, tmp)
  ssf = np.frombuffer(tmp, dtype=np.double, count=npoints * ncol)
  return ssf.reshape((npoints, ncol))
