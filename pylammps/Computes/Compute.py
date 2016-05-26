"""
Compute class
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
#TODO: Take care of the looks of these plots

class Compute(object):
  """
  Abstract compute class. It will never be used, but is parent of all
  the different computes.
  """
  def __init__(self):
    """
    Constructor. Not clear what to do here
    """
    self.value = 0
    self.idx = 0
    self.header = []

  def compute(self, system):
    """
    Compute routine
    """
    pass

  def tally(self, value):
    """
    Tally new compute with the previous ones. Mostly because not all
    of the computes have the same structure, so the "average" is not
    standard. By default we do the usual average.
    """
    self.idx += 1
    self.value *= (self.idx - 1)/self.idx
    self.value += value/self.idx

  def zero(self):
    """
    Zero out current tallies.
    """
    self.value = 0
    self.idx = 0

  def log(self, filename):
    """
    Logging routine. By default we just write self.value to filename,
    with self.header
    """
    np.savetxt(filename, self.value, header='; '.join(self.header))

  def plot(self, filename):
    """
    Plotting routine. By default we plot every column [1:] as a
    function of column 0, setting labels and axis names with
    self.header and save it to filename.
    """
    fig, axis = pl.subplots()
    for i, vec in enumerate(self.value.T[1:]):
      axis.plot(self.value[:, 0], vec, label=self.header[i])
    axis.set_xlabel(self.header[0])
    fig.savefig('{0}.pdf'.format(filename))
    pl.close()
