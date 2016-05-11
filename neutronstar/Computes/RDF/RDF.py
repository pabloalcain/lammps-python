"""
Compute class
"""

from neutronstar.Computes import Compute
import ctypes as ct
import numpy as np
import os
#TODO: Take care of the looks of these plots
_DIRNAME = os.path.dirname(__file__)
librdf = ct.CDLL(os.path.join(_DIRNAME, 'librdf.so'))
rdf_c = librdf.rdf

class RDF(Compute):
  """
  Radial Distribution Function calculation, order N^2 computing *all
  neighbors*.
  """
  def __init__(self, nbins, pairs, pbc=True):
    """
    Constructor.

    Parameters
    ----------

    nbins : int
        Number of bins to consider

    pairs: iterable
        List of pairs to consider. Each element of the list is a tuple
        of size 2, with the first element being the types considered
        and the second one the other types. For example, if we want to
        calculate the RDF on a system with 2 types of particles of all
        types vs all types and type 1 vs type 2, the pairs should be:

        pairs = [((1, 2), (1, 2)), ((1,), (2,))]

        A keyword for "all types" is type 0

    pbc: boolean, optional
        Whether or not to use periodic boundary conditions. Default is
        True.

    Note
    ----

    So far we cannot put more than one key inside each tuple, since it
    would probably either some serious refactoring or messy code. 'All
    particles' can still be chosen by keyword value 0
    """
    self.pairs = pairs
    self.nbins = nbins
    self.pbc = pbc
    #self.rdf = np.zeros((nbins, len(pairs) + 1))
    super(RDF, self).__init__()
    self.header = []
    for p in pairs:
      h = []
      for t in p:
        h.append(','.join(map(str, t)))
      self.header.append('-'.join(map(str, h)))

  def compute(self, system):
    """Calculate RDF.

    Parameters
    ----------

    system : System
        System on which we calculate the Radial Distribution Function

    Returns
    -------

    rdf : numpy array
        An array with the information of the compute. The first column
        is 'r' and the rest are the RDF calculated for the pair list.
    """
    #TODO: Add some warning in case the indices repeat
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
    tmp = (ct.c_double * (self.nbins * ncol))()
    x_p = system.x.ctypes.data_as(ct.c_void_p)
    t_p = system.t.ctypes.data_as(ct.c_void_p)
    pair_p = pair_ar.ctypes.data_as(ct.c_void_p)
    rdf_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int, ct.c_int,
                      ct.c_double, ct.c_void_p, ct.c_int, ct.c_bool,
                      ct.c_void_p]
    rdf_c(x_p, t_p, natoms, self.nbins, size, pair_p, npairs,
          self.pbc, tmp)
    rdf = np.frombuffer(tmp, dtype=np.double, count=self.nbins * ncol)
    return rdf.reshape((self.nbins, ncol))
