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
import matplotlib.pyplot as plt
import numpy as np

_tags = {'breadth': 'Breadth [fm]',
         'del_lambda': r'Fitting $\Delta\lambda$ [fm]',
         'energy': 'Energy [MeV]',
         'pressure': r'Pressure [MeV/fm$^3$]',
         'S_absorption': r'Pasta absorption',
         'size_avg': 'Average size of clusters',
         'k_absorption': r'Pasta wavenumber [fm$^{-1}$]',
         'surface': r'Surface [fm$^2$]',
         'height': r'Fitting height',
         'volume': r'Volume [fm$^3$]',
         'size_std': r'$\sigma$ cluster distribution',
         'potential': r'Potential energy [MeV]',
         'temperature': r'Temperature [MeV]',
         'kinetic': r'Kinetic energy [MeV]',
         'del_height': r'Fitting $\Delta$ H',
         'euler': r'Euler characteristic',
         'lambda': r'Fitting wavelength [fm]',}
    


def mste((occ, frac), path='./'):
    indices = (occ != 0.0)
    mass = np.arange(len(occ)) # mass is just the indices of occ and frac
    frac = frac[indices]
    mass = mass[indices]
    occ = occ[indices] # Filter out
    np.savetxt(path + 'cluster.dat', np.transpose([mass, occ, frac]),
               header='size, number, frac', fmt='%i, %f, %f')

    # Coarse-grain the histogram
    nbins = 20
    bins = np.logspace(0, np.log10(mass[-1]), nbins)
    a, xx = np.histogram(mass, weights = occ/50.0, bins = bins)

    # Set the figure
    fig, ax1 = plt.subplots()
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
    plt.close(fig)

def rdf(rdf, path='./'):
    
    np.savetxt(path + 'rdf.dat', rdf,
               header=('r, a-a, ia-a, n-n, in-n, '
                       'n-p, in-p, p-p, ip-p'),
               fmt=', '.join(['%f']*9))

    labels = ('all', 'n-n', 'n-p', 'p-p')
    fig, ax = plt.subplots()
    ax.set_xlabel('Distance [fm]')
    ax.set_ylabel('RDF')
    R = rdf[:, 0]
    ax.set_xlim(0, R[-1])
    for i in xrange(4):
        pair = labels[i]
        idx = i * 2 + 1
        G = rdf[:, idx]
        fig_each, ax_each = plt.subplots()
        ax_each.set_xlabel('Distance [fm]')
        ax_each.set_xlim(0, R[-1])
        for axis in (ax, ax_each):
            axis.plot(R, G, label=pair)
        ax_each.set_ylabel('RDF ({0})'.format(pair))
        fig_each.tight_layout()
        fig_each.savefig(path + '/rdf_{0}.pdf'.format(pair))
        plt.close(fig_each)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path + '/rdf.pdf')
    plt.close(fig)

def ssf(ssf, path='./'):

    np.savetxt(path + 'ssf.dat', ssf,
               header=('r, a-a, n-n, n-p, p-p'),
               fmt=', '.join(['%f']*5))

    labels = ('all', 'n-n', 'n-p', 'p-p')
    fig, ax = plt.subplots()
    ax.set_xlabel(r'Wave number [fm$^{-1}$]')
    ax.set_ylabel('SSF')
    Q = ssf[:, 0]
    qmax = 7.0
    ax.set_xlim(0, qmax)
    for i in xrange(4):
        idx = i + 1
        pair = labels[i]
        S = ssf[:, idx]
        fig_each, ax_each = plt.subplots()
        ax_each.set_xlabel(r'Wave number [fm$^{-1}$]')
        ax_each.set_xlim(0, qmax)
        for axis in (ax, ax_each):
            axis.plot(Q, S, label=pair)
        ax_each.set_ylabel('SSF ({0})'.format(pair))
        fig_each.tight_layout()
        fig_each.savefig(path + '/ssf_{0}.pdf'.format(pair))
        plt.close(fig_each)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path + '/ssf.pdf')
    plt.close(fig)

def thermo(collective, path='./'):

    data = np.vstack(collective.values()).T
    header = ', '.join(collective.keys())
    np.savetxt(path + 'thermo.dat', data, header=header, fmt='%f')

    for i, key in enumerate(collective.keys()):
        fig, ax = plt.subplots()
        ax.set_xlabel('Timestep')
        ax.set_ylabel(_tags[key])
        ax.plot(collective[key])
        fig.tight_layout()
        fig.savefig(path + key + '.pdf')
        plt.close(fig)
