"""
This automates many plots usually done when post-processing
"""

import matplotlib
params = {'axes.labelsize': 12,
          'font.size': 15,
          'legend.fontsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'text.usetex': False,
          'lines.linewidth': 2.0,
          'lines.markersize': 10.0,
          'legend.numpoints': 1,}

def _format(axis):
  """
  Format axis to make them pretty
  """
  axis.spines['top'].set_visible(False)
  axis.spines['right'].set_visible(False)
  axis.spines['left'].set_visible(False)
  axis.get_xaxis().tick_bottom()
  axis.get_yaxis().tick_left()
  axis.tick_params(axis='x', direction='out')
  axis.tick_params(axis='y', length=0)
  for spine in axis.spines.values():
    spine.set_position(('outward', 5))
  axis.grid(axis='y', color='0.9', linestyle='-', linewidth=1)
  axis.set_axisbelow(True)




matplotlib.rcParams.update(params)
from itertools import cycle

import numpy as np
import matplotlib.pyplot as plt
from extract import Extraction
import analysis as A

COLOR = ('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854',
         '#ffd92f', '#e5c494', '#b3b3b3')
#COLOR = ('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
#         '#e6ab02', '#a6761d', '#666666')


FMTL = ['^', 's', 'D', 'v', 'o',]
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

    txmin = min(xval)
    txmax = max(xval)
    tymin = min(yval - yerr/2)
    tymax = max(yval + yerr/2)
    if len(axis.get_lines()) == 0:
      axis.set_xlim((txmin, txmax))
      axis.set_ylim((tymin, tymax))
    else:
      xmin, xmax = axis.get_xlim()
      axis.set_xlim((min(xmin, txmin), max(xmax, txmax)))
      ymin, ymax = axis.get_ylim()
      axis.set_ylim((min(ymin, tymin), max(ymax, tymax)))

    axis.set_xlabel(LABELS[var])
    axis.set_ylabel(LABELS[mag])

    axis.plot(xval, yval, **kwargs)
    axis.fill_between(xval, yval - yerr/2, yval + yerr/2, alpha=0.1,
                      color=kwargs['color'])
    return

  def multiplot(self, run, mag, parameters):
    """
    From a parameters dict that has two iterables, select one as
    'running variable' and plot for different values of it
    """
    _param = parameters.copy()
    fig, axis = plt.subplots(1, 1)
    for fmt, color, val in zip(FMTL, COLOR, parameters[run]):
      _param[run] = val
      lbl = '{0} = {1}'.format(run, val)
      self.onedim(axis, mag, _param, marker=fmt, label=lbl, color=color)
    legend = axis.legend(loc='upper left')
    _format(axis)
    frame = legend.get_frame()

    frame.set_facecolor('1.0')
    frame.set_edgecolor('1.0')
    fig.canvas.draw()
    fig.tight_layout()
    return fig, axis

  def cluster(self, parameters, load=True, idx=0):
    """Plot cluster of parameters dict. If load, extract cluster from .dat
    file, else recompute. If we are going to recompute, choose idx
    of the lammpstrj file to recompute.
    """
    if load:
      cl = self.extraction.array('cluster', parameters)
    else:
      index = self.extraction.mste(parameters, idx=idx)
      t = self.extraction.particle('type', parameters, idx=idx)
      cl = A.reduce_mste(index, t)
    mass = cl[:, 0]
    occ = cl[:, 1]
    frac = cl[:, 2]
    nbins = 20
    bins = np.logspace(0, np.log10(mass[-1]), nbins)

    a, xx = np.histogram(mass, weights=occ, bins=bins)

    # "Averaging" the proton fraction inside each histogram
    frac2, _ = np.histogram(mass, weights=occ*frac, bins=bins)
    mask = (a != 0)
    frac2 = frac2[mask]/a[mask]
    xm = (xx[:-1] + xx[1:])/2
    xm = xm[mask]
    # Set the figure
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Cluster size')
    ax1.set_xscale('log')
    ax2 = ax1.twinx()
    # Plot frequency
    ax1.bar(xx[:-1], a, np.diff(xx), alpha=0.5, label='Frequency', color=COLOR[0])
    ax1.set_yscale('log')
    ax1.set_ylabel('Frequency', color=COLOR[0])
    h1, l1 = ax1.get_legend_handles_labels()
    for t1 in ax1.get_yticklabels():
      t1.set_color(COLOR[0])
    # Plot proton fraction
    ax2.plot(xm, frac2, 'o--', label='Proton fraction', color=COLOR[1])
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Proton fraction', color=COLOR[1])
    h2, l2 = ax2.get_legend_handles_labels()
    for t2 in ax2.get_yticklabels():
      t2.set_color(COLOR[1])
    # Add legend save and close fig
    ax1.legend(h1+h2, l1+l2)
    fig.canvas.draw()
    fig.tight_layout()
    return fig, ax1

  def rdf(self, parameters, pair='n-n', load=True, idx=0):
    """Plot rdf of parameters dict. If load, extract cluster from .dat
    file, else recompute. If we are going to recompute, choose idx
    of the lammpstrj file to recompute.
    """
    _d = {'all': 1, 'p-p': 2, 'n-p': 3, 'n-n': 4}

    if load:
      rdf = self.extraction.array('rdf', parameters)
    else:
      rdf = self.extraction.rdf(parameters, idx=idx)
    fig, axis = plt.subplots(1, 1)
    i = _d[pair]
    axis.plot(rdf[:, 0], rdf[:, i], color=COLOR[0])
    axis.set_xlabel('Distance [fm]')
    axis.set_ylabel('RDF')
    fig.canvas.draw()
    fig.tight_layout()
    return fig, axis


def main():
  """Main function with a use case, where we plot the EoS"""
  plot = Plotter('/home/pablo/artemis/eos/data')
  parameters = {'x': (0.5, 0.4, 0.3, 0.2, 0.1),
                'd': np.linspace(0.02, 0.16, 15),
                'T': 0.5, 'V': 'medium'}
  run = 'x'
  mag = 'energy'
  plot.multiplot(run, mag, parameters)
  parameters = {'x': 0.4, 'd': 0.04, 'T': 1.5, 'V': 'medium'}
  plot.cluster(parameters, load=False, idx=1)
  plot.rdf(parameters, load=False)
  parameters = {'x': 0.5,
                'd': np.linspace(0.02, 0.16, 15),
                'T': (0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0),
                'V': 'medium'}
  run = 'T'
  mag = 'lambda'
  plot.multiplot(run, mag, parameters)

  plt.show()

if __name__ == '__main__':
  main()
