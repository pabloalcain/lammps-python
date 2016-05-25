"""
Compute class
"""

from neutronstar.Computes import Compute
from analysis import rdf

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

    pbc : boolean, optional
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
    #this should not be needed
    #self.rdf = np.zeros((nbins, len(pairs) + 1))
    super(RDF, self).__init__()
    self.header = []
    for p in pairs:
      h = []
      try:
        for t in p:
          h.append(','.join(map(str, t)))
      except TypeError:
        h.append(str(p))
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
    val = rdf(system.x, system.t, system['size'], self.pairs,
              self.nbins, pbc=self.pbc)
    return val
