"""
Analyzer class: the Analyzer of LAMMPS Systems
"""

class Analyzer(object):
  def __init__(self, computes):
    """
    Constructor: Instantiate an analyzer with a compute list

    Parameters
    ----------

    computes : dict
        A dict of computes that will be used to analyze. The key is
        the identifier and the value is a compute instance
    """
    pass

  @property
  def computes(self):
    """
    Computes for the analyzer
    """
    pass

  @computes.setter
  def computes(self, value):
    """
    Set computes for the analyzer
    """
    pass

  def analyze(self, system):
    """
    Main method, that performs all the analysis on the system given

    Parameters
    ----------

    system : System instance
        The system that we will use to calculate the computes
    """
    pass
