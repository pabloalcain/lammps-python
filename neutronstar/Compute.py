"""
Compute class
"""

class Compute(object):
  """
  Abstract compute class. It will never be used, but is parent of all
  the different computes.
  """
  def __init__(self):
    """
    Constructor. Not clear what to do here
    """
    pass

  def compute(self):
    """
    Compute routine
    """
    pass

  def tally(self):
    """
    Tally new compute with the previous ones. Mostly because not all
    of the computes have the same structure, so the "average" is not
    standard.
    """
    pass

  def zero(self):
    """
    Zero out current tallies.
    """
    pass

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
    """
    pass

  def compute(self, system):
    """
    Calculate RDF.

    Parameters
    ----------

    system : System
        System on which we calculate the Radial Distribution Function

    Returns
    -------

    rdf : dict
        A dictionary with the information of the compute. It has as
        keys 'r' and the list of pairs in which the RDF was
        calculated.
    """
    pass

  def tally(self, system):
    """
    """

class MST(Compute):
  """
  MST/MSTE calculation.
  """
  def __init__(self, energy=True):
    """
    Constructor.

    Parameters
    ----------

    energy : boolean
        Whether or not to use energy considerations (MSTE) in the
        cluster recognition. Default is True.
    """
    pass

  def compute(self, system, pairs):
    """
    Calculate MST.

    Parameters
    ----------

    system : System
        System on which we calculate the Minimum Spanning Tree

    Returns
    -------

    mste : list
        A list of tuples, with format (mass, x, infinite): mass is the
        total mass of the system, x is the proton fraction and
        infinite is boolean on whether the fragment is infinite or
        not.
    """
    pass

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
    """
    pass

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
        A dictionary with the information of the compute. It has as
        keys 'k' and the list of pairs in which the structure factor
        was calculated.
    """
    pass
