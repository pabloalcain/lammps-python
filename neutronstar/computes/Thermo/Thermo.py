"""
Thermodynamic Compute class
"""

from neutronstar.computes.Compute import Compute

class Thermo(Compute):
  """
  MST/MSTE calculation.
  """
  def __init__(self):
    """
    Constructor.

    Parameters
    ----------

    energy : boolean
        Whether or not to use energy considerations (MSTE) in the
        cluster recognition. Default is True.
    """
    super(Thermo, self).__init__()
    self.header = ['Temperature', 'Kinetic', 'Potential', 'Total',
                   'Pressure']
    self.value = []

  def compute(self, system):
    """Calculate MSTE.

    Parameters
    ----------

    system : System
        System on which we calculate the Minimum Spanning Tree

    Returns
    -------

    value, (mst, inf) : numpy array, numpy array, list
        value is the [mass, occupancy, fraction] histogram
        mst is the array of indices to which each particle belongs.
        inf is the list of infinite clusters.

    Notes
    -----

    So far, when energy is `True`, this only works for systems with
    `medium` potential.
    """
    natoms = system['N']
    temperature = system.lmp.extract_compute("thermo_temp", 0, 0)
    potential = system.lmp.extract_compute("thermo_pe", 0, 0)/natoms
    pressure = system.lmp.extract_compute("thermo_press", 0, 0)
    kinetic = 3.0/2.0 * temperature
    return temperature, kinetic, potential, potential + kinetic, pressure

  def tally(self, value):
    """
    We need to override the parent tally, since this compute does
    not average trivially. We create the histogram from the values.
    """
    self.idx += 1
    self.value.append(value)

  def zero(self):
    """
    The value has to be an empty list
    """
    super(Thermo, self).zero()
    self.value = []
