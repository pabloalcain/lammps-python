from scipy.stats.distributions import t
import fit as F
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
    p = self.set_path(V, l, x, N, d, T)
    fp = open(p + "/thermo.dat")
    l = fp.readline()[2:-1]
    fp.close()
    
    try:
      idx = l.split(', ').index(mag)
    except ValueError:
      if mag in ['lambda', 'height']:
        idx = -1
      else:
        raise ValueError

    if idx != -1:
      A = np.loadtxt(p +'/thermo.dat')
      c = np.mean(A, axis=0)[idx]
      s = np.std(A, axis=0)[idx]
    else:
      sm_popt, sm_pcov = F.wavelength(p + '/rdf.dat')
      tval = t.ppf(1.0-0.31/2, 3)
      if mag == "lambda":
        s = sm_pcov[1][1]**0.5 * tval
        c = sm_popt[1]
      if mag == "height":
        s = sm_pcov[0][0]**0.5 * tval
        c = sm_popt[0]
    return c, s

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
