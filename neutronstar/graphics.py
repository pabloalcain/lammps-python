"""
Create plots and raw data for different magnitudes usual in
neutronstars
"""

import matplotlib
font = {'size' : 18}
matplotlib.rc('font', **font)
lines = {'linewidth' : 2.0, 'markersize': 10.0}
matplotlib.rc('lines', **lines)
legend = {'numpoints': 1}
matplotlib.rc('legend', **legend)
matplotlib.use('Agg')
import pylab as pl
import numpy as np

def flush_mste(occ, frac, path='./'):
    indices = (occ != 0.0)
    mass = np.arange(len(occ)) # mass is just the indices of occ and frac
    frac = frac[indices]
    mass = mass[indices]
    occ = occ[indices] # Filter out
    np.savetxt(path + 'cluster.dat', np.transpose([mass, occ, frac]),
               header='size, number, frac', fmt='%i, %f, %f')

    # Coarse-graom the histogram
    nbins = 20
    bins = np.logspace(0, np.log10(mass[-1]), nbins)
    a, xx = np.histogram(mass, weights = occ/50.0, bins = bins)

    
    fig, ax1 = pl.subplots()
    ax1.set_xlabel('Cluster size')
    ax1.set_xscale('log')
    ax2 = ax1.twinx()
    # Plot frequency
    ax1.bar(xx[:-1], a, np.diff(xx), alpha=0.5, label='Frequency')
    ax1.set_yscale('log')
    ax1.set_ylabel('Frequency', color='b')
    h1, l1 = ax1.get_legend_handles_labels()
    for t1 in ax1.get_yticklabels(): t1.set_color('b')
    # Plot proton fraction
    ax2.plot(mass, frac, 'r^--', label='Proton fraction')
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Proton fraction', color='r')
    h2, l2 = ax2.get_legend_handles_labels()
    for t2 in ax2.get_yticklabels(): t2.set_color('r')
    # Add legend save and close fig
    ax1.legend(h1+h2, l1+l2)
    fig.tight_layout()
    fig.savefig(path + 'cluster.pdf')
    pl.close(fig)
