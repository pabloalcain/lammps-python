import matplotlib
from matplotlib import pyplot as plt
import numpy as np
params = {'axes.labelsize': 12,
          'font.size': 15,
          'legend.fontsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'text.usetex': False,
          'lines.linewidth': 2.0,
          'lines.markersize': 10.0,
          'legend.numpoints': 1,
          'axes.color_cycle': ('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3',
                               '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3')}

#COLOR = ('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
#         '#e6ab02', '#a6761d', '#666666')


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

import postprocess

def mste(value, energy=False):
  fig, ax1 = plt.subplots()
  ax1.set_xlabel('Cluster size')
  ax1.set_xscale('log')
  ax2 = ax1.twinx()
  if energy:
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("axes", 1.2))
  mass = value[1:, 0]
  occ = value[1:, 1]
  frac = value[1:, 2]
  if energy: clus_en = value[1:, 3]

  nbins = 30
  bins = np.logspace(0, np.log10(mass[-1]), nbins)

  a, xx = np.histogram(mass, weights=occ, bins=bins)

  # "Averaging" the proton fraction inside each histogram
  frac2, _ = np.histogram(mass, weights=occ*frac, bins=bins)
  mask = (a != 0)
  frac2 = frac2[mask]/a[mask]
  xm = (xx[:-1] + xx[1:])/2
  xm = xm[mask]
  if energy:
  # "Averaging" the proton fraction inside each histogram
    clus_en2, _ = np.histogram(mass, weights=occ*clus_en, bins=bins)
    mask = (a != 0)
    clus_en2 = clus_en2[mask]/a[mask]

  ax1.plot(mass[occ > 0], occ[occ > 0], '.-', label='Frequency')
  ax1.set_yscale('log')
  ax1.set_ylabel('Frequency')
  hand1, lab1 = ax1.get_legend_handles_labels()
  # Plot proton fraction
  color_cycle = ax1._get_lines.color_cycle
  ax2.plot(xm, frac2, 'o--',
           label='Proton fraction', color=next(color_cycle))
  ax2.set_ylim(0, 1)
  ax2.set_ylabel('Proton fraction')
  hand2, lab2 = ax2.get_legend_handles_labels()
  if energy:
    ax3.plot(xm, clus_en2, 'o--',
             label='Energy', color=next(color_cycle))
    ax3.set_ylabel('Energy')
    hand3, lab3 = ax3.get_legend_handles_labels()
  # Add legend save and close fig
  if energy:
    ax1.legend(hand1+hand2+hand3, lab1+lab2+lab3)
  else:
    ax1.legend(hand1+hand2, lab1+lab2)
  fig.tight_layout()
  _format(ax1)
  return fig, ax1, ax2
