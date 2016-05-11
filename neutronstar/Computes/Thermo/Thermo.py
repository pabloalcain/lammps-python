"""
Thermodynamic Compute class
"""

from neutronstar.Computes.Compute import Compute
import numpy as np
import pylab as pl

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
    """Calculate Thermo.

    Parameters
    ----------

    system : System
        System on which we calculate the Minimum Spanning Tree

    Returns
    -------

    temperature, kinetic, potential, total, pressure: 5*float
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

  def plot(self, filename):
    """
    We need to override the parent plot method
    """
    work = np.array(self.value).T
    for i, measure in enumerate(work):
      fig, axis = pl.subplots()
      axis.plot(measure)
      axis.set_ylabel(self.header[i])
      fig.savefig('{0}-{1}.pdf'.format(filename, self.header[i]))
      pl.close()
