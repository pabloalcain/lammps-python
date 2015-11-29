"""
This automates many plots usually done when post-processing
"""

import matplotlib
_FONT = {'size'   : 18}
matplotlib.rc('font', **_FONT)
_LINES = {'linewidth' : 2.0, 'markersize': 10.0}
matplotlib.rc('lines', **_LINES)
_LEGEND = {'numpoints': 1}
matplotlib.rc('legend', **_LEGEND)
import numpy as np
import matplotlib.pyplot as plt
from extract import Extraction

FMTL = ['^-', 's-', 'D-', 'v-', 'o-',]
LABELS = {'energy': r'Energy [MeV]',
          'lambda': r'$\lambda$ [fm]',
          'S_absorption': 'Absorption peak',
          'pressure': r'Pressure [MeV/fm$^3$]',
          'euler': r'Euler Number',
          'breadth': r'Mean Breadth [fm]',
          'surface': r'Mean Surface [fm$^2$]',
          'volume': r'Mean Volume [fm$^3$]',
          'd': r'Density [fm$^{-3}$]',
          'T': r'Temperature [MeV]',}

SRC = 'source_fig/'
DEST = 'fig/'

class Plotter(object):
  """
  We get a path and instantiate Extraction. From then on, we plot as
  abstractly as possible.

  # TODO: Add plots for arrays (like mste, rdf, etc)
  """
  def __init__(self, path='.'):
    self.extraction = Extraction(path)

  def onedim(self, axis, mag, parameters, **kwargs):
    """
    This function plots in axis a magnitude mag as a function of
    whichever parameter that is an interable
    """
    _param = parameters.copy()
    for key, value in parameters.iteritems():
      try:
        len(value)
        if not isinstance(value, str):
          var = key
          xval = value
      except TypeError:
        pass

    yval = np.zeros_like(xval)
    yerr = np.zeros_like(xval)
    for i, x in enumerate(xval):
      _param[var] = x
      ext = self.extraction.scalar(mag, _param)
      yval[i] = ext[0]
      yerr[i] = ext[1]
    axis.errorbar(xval, yval, yerr=yerr, **kwargs)
    axis.set_xlabel(LABELS[var])
    axis.set_ylabel(LABELS[mag])
    return

  def multiplot(self, run, mag, parameters):
    """
    From a parameters dict that has two iterables, select one as
    'running variable' and plot for different values of it
    """
    _param = parameters.copy()
    fig, axis = plt.subplots(1, 1)
    for fmt, val in zip(FMTL, parameters[run]):
      _param[run] = val
      lbl = '{0} = {1}'.format(run, val)
      self.onedim(axis, mag, _param, fmt=fmt, label=lbl)
    axis.legend()
    fig.canvas.draw()
    fig.tight_layout()
    return fig, axis

def main():
  """Main function with a use case, where we plot the EoS"""
  plot = Plotter('/home/pablo/artemis/eos/data')
  parameters = {'x': (0.1, 0.2, 0.3, 0.4, 0.5),
                'd': np.linspace(0.02, 0.16, 15),
                'T': 0.5, 'V': 'medium'}
  run = 'x'
  mag = 'eos'
  plot.multiplot(run, mag, parameters)
  plt.show()

if __name__ == '__main__':
  main()
