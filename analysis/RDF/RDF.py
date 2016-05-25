"""
RDF calculation
"""

import ctypes as ct
import numpy as np
import os

_DIRNAME = os.path.dirname(__file__)
librdf = ct.CDLL(os.path.join(_DIRNAME, 'librdf.so'))
rdf_c = librdf.rdf

def rdf(x, t, size, pairs, nbins, pbc=True):
  """Calculate RDF.

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

  nbins : int
      Number of bins to consider

  pbc : boolean, optional
      Whether or not to use periodic boundary conditions. Default is
      True.

  Returns
  -------

  rdf : numpy array
      An array with the information of the compute. The first column
      is 'r' and the rest are the RDF calculated for the pair list.
  """
  #TODO: Add some warning in case the indices repeat
  natoms = np.shape(x)[0]
  npairs = len(pairs)
  ncol = npairs + 1
  pair_ar = np.zeros(2*npairs)
  i = 0
  for p in pairs:
    pair_ar[i] = p[0][0]
    pair_ar[i+1] = p[1][0]
    i += 2
  tmp = (ct.c_double * (nbins * ncol))()
  x_p = x.ctypes.data_as(ct.c_void_p)
  t_p = t.ctypes.data_as(ct.c_void_p)
  pair_p = pair_ar.ctypes.data_as(ct.c_void_p)
  rdf_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int, ct.c_int,
                    ct.c_double, ct.c_void_p, ct.c_int, ct.c_bool,
                    ct.c_void_p]
  rdf_c(x_p, t_p, natoms, nbins, size, pair_p, npairs, pbc, tmp)
  value = np.frombuffer(tmp, dtype=np.double, count=nbins * ncol)
  return value.reshape((nbins, ncol))
