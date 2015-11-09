#from scipy.stats.distributions import t
import numpy as np
#import fit as F
import os

class Plotter(object):
  def __init__(self, path = '.'):
    self.path = path

  def check(self, p):
    if len(os.listdir(p)) > 1:
      helper = ("Multiple options for directory {0}: "
                "{1}. Please specify a value")
      raise AttributeError(helper.format(p,
                                         os.listdir(p)))
    return os.listdir(p)[0]

  def set_path(self, V, l, x, N, d, T):
    #Here we believe on the dir structure:
    #{V}/l{l}/x{x}/N{N}/d{d}/T{T}

    p = self.path

    if V == None:
      direct = self.check(p)
    else:
      direct = str(V)
    p += "/" + direct

    if l == None:
      direct = self.check(p)
    else:
      direct = "l" + str(l)
    p += "/" + direct

    if x == None:
      direct = self.check(p)
    else:
      direct = "x" + str(x)
    p += "/" + direct

    if N == None:
      direct = self.check(p)
    else:
      direct = "N" + str(N)
    p += "/" + direct

    if d == None:
      direct = self.check(p)
    else:
      direct = "d" + str(d)
    p += "/" + direct

    if T == None:
      direct = self.check(p)
    else:
      direct = "T" + str(T)
    p += "/" + direct
    return p


  def extract(self, mag, V = None, 
              l = None, x = None, 
              N = None, d = None, 
              T = None):
    """
    It doesn't support anymore height, that came from the
    fit method. Now lambda makes the conversion from k to lambda.
    """
    p = self.set_path(V, l, x, N, d, T)
    try:
      fp = open(p + "/thermo.dat")
    except IOError:
      return 0, 0
    l = fp.readline()[2:-1]
    #    l = "breadth, del_lambda, S_absorption, pressure, size_avg, k_absorption, surface, height, volume, size_std, del_height, potential, kinetic, lambda, energy, euler, temperature"
    fp.close()

    if mag == 'lambda':
      idx = l.split(', ').index('k_absorption')
    else:
      idx = l.split(', ').index(mag)
    A = np.loadtxt(p +'/thermo.dat')
    c = np.mean(A, axis=0)[idx]
    s = np.std(A, axis=0)[idx]
    if mag == 'lambda': 
      s = 2*np.pi/(c*c) * s
      c = 2*np.pi/c
    return c, s

  def extract_array(self, mag, V = None, 
                    l = None, x = None, 
                    N = None, d = None, 
                    T = None):
    p = self.set_path(V, l, x, N, d, T)
    try:
      A = np.loadtxt('{0}/{1}.dat'.format(p, mag))
    except ValueError:
      A = np.loadtxt('{0}/{1}.dat'.format(p, mag), delimiter = ',')
    return A
  
  

if __name__ == "__main__":
  import numpy as np
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as pl
  p = Plotter('/home/pablo/lammps-neutron-stars-input/examples/data_explor')
  mag = 'euler'
  
  x_l = np.linspace(0.4,0.5,6)
  d_l = np.linspace(0.01,0.08,8)
  T_l = np.linspace(0.5, 2.0, 31)
  for x in x_l:
    for d in d_l:
      med = []
      stiff = []
      for T in T_l:
        out = p.extract(mag, V="medium", x=x, d=d, T=T)
        med.append(out)
        out = p.extract(mag, V="stiff", x=x, d=d, T=T)
        stiff.append(out)
      med = np.array(med)
      stiff = np.array(stiff)

      pl.figure()
      ax = pl.subplot(1,1,1)
      ax.errorbar(T_l, med[:, 0], marker="o",
                  yerr=med[:, 1], label="medium")
      ax.errorbar(T_l, stiff[:, 0], marker="o",
                  yerr=stiff[:, 1], label="stiff")
      ax.legend(numpoints=1)
      pl.xlabel('Temperature [MeV]')
      pl.ylabel(mag)
      pl.savefig('{mag}-x{x}-d{d}.pdf'.format(mag=mag,
                                              x=x,
                                              d=d))
