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
    for key in self.computes:
      self.computes[key].compute(system)

  def update(self, system):
    """
    Perform the computes and update all the tallies

    Parameters
    ----------

    system : System instance
        The system that we will update
    """
    for key in self.computes:
      value = self.computes[key].compute(system)
      self.computes[key].tally(value)

  def zero(self):
    """
    Zero out all the computes
    """
    for key in self.computes:
      self.computes[key].zero()

  def log(self, path):
    """
    Log all the computes

    Parameters
    ----------

    path : str
        Path where to write the log

    """
    for key in self.computes:
      self.computes[key].log('{0}/{1}.dat'.format(path, key))

  def plot(self, path):
    """
    Plot all the computes

    Parameters
    ----------

    path : str
        Path where to save the plots

    """
    for key in self.computes:
      self.computes[key].plot('{0}/{1}.pdf'.format(path, key))
