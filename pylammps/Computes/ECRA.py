"""
ECRA Compute class
"""

from pylammps.Computes import Compute
from analysis import ecra
import pylab as pl

class ECRA(Compute):
  """
  ECRA calculation.
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
    super(ECRA, self).__init__()
    self.header = ['mass', 'occupancy', 'fraction']

  def compute(self, system):
    """Calculate ECRA.

    Parameters
    ----------

    system : System
        System on which we calculate the Minimum Spanning Tree

    Returns
    -------

    value : numpy array
        value is the [mass, occupancy, fraction] histogram
    """
    val = ecra(system.x, system.v, system.t, system.box,
               system['expansion'])
    return val[0]

  def tally(self, value):
    """
    We need to override the parent tally, since this compute does
    not average trivially. We create the histogram from the values.
    """
    self.idx += 1
    #If we are in the first frame of the calculation, just replace by
    #the value.
    try:
      new_occ = (self.idx - 1) * self.value[:, 0] + value[:, 0]
      self.value[:, 0] = new_occ / self.idx
      self.value[:, 1] *= (self.idx - 1) * self.value[:, 0]
      self.value[:, 1] += value[:, 1] * value[:, 0]
      self.value[:, 1] /= new_occ
    except TypeError:
      if self.idx == 1:
        self.value = value
      else: raise TypeError

  def plot(self, filename):
    """
    Plotting routine. We need to override since in this case we don't
    want both plots to be on the same y-axis
    """
    fig, ax1 = pl.subplots()
    ax1.set_xlabel('Cluster size')
    ax1.set_xscale('log')
    ax2 = ax1.twinx()
    mass = self.value[1:, 0]
    occ = self.value[1:, 1]
    frac = self.value[1:, 2]
    ax1.plot(mass[occ > 0], occ[occ > 0], '.-', label='Frequency')
    ax1.set_yscale('log')
    ax1.set_ylabel('Frequency')
    hand1, lab1 = ax1.get_legend_handles_labels()
    # Plot proton fraction
    ax2.plot(mass[occ > 0], frac[occ > 0], 'o--', label='Proton fraction')
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Proton fraction')
    hand2, lab2 = ax2.get_legend_handles_labels()
    # Add legend save and close fig
    ax1.legend(hand1+hand2, lab1+lab2)
    fig.tight_layout()
    fig.savefig('{0}.pdf'.format(filename))
    pl.close()
