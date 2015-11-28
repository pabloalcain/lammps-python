"""
This is simply a wrapper for the library libanalysis.so
"""

import ctypes as ct
import numpy as np
from lammps import lammps
from scipy.optimize import curve_fit
from scipy.stats.distributions import t

analysis = ct.cdll.LoadLibrary("libanalysis.so")

def rdf(lmp, nbins, npairs):
  """
  Wrapper for the c++ libanalysis. Its main disadvantage is that we
  should know beforehand the number of pairs to consider.  Since
  this is meant to Neutron Stars, we document here behavior when
  npair is set to 4:

  If we label the pairs with index k, so that
  k = 1 => all v all
  k = 2 =>  1  v  1
  k = 3 =>  1  v  2
  k = 4 =>  2  v  2

  Returns a 2D array, with this format.

  1st column: Distance
  k + 1: g(r) between the pair labeled by k
  """
  x = lmp.gather_atoms('x', 1, 3)
  typ = lmp.gather_atoms('type', 0, 1)
  natoms = lmp.get_natoms()
  dx = lmp.extract_global('boxxhi', 1) - lmp.extract_global('boxxlo', 1)
  dy = lmp.extract_global('boxyhi', 1) - lmp.extract_global('boxylo', 1)
  dz = lmp.extract_global('boxzhi', 1) - lmp.extract_global('boxzlo', 1)
  if not dx == dy == dz:
    raise ValueError('Cannot compute for non-cubic boxes! exiting')
  size = dx
  ncol = npairs + 1
  tmp = (ct.c_double * (nbins * ncol))()
  analysis.rdf.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int, ct.c_int, ct.c_double, ct.c_void_p]
  analysis.rdf(x, typ, natoms, nbins, size, tmp)
  r = np.frombuffer(tmp, dtype=np.float, count=nbins * ncol)
  return np.reshape(r, (nbins, ncol))

def structure(gr, density, npairs):
  """
  Calculate structure factor given the radial distribution function.
  Returns an array structured like rdf, with only one column per
  pair.

  Also return the first peak of the neutron-neutron scattering.
  (peak between ~8 and ~30 fm wavelength)
  """
  _d = density
  r = gr[:, 0]
  # Assume evenly spaced
  dr = r[1] - r[0]
  # How many points do I need to add to get a good resolution?
  dlda = 0.1
  lda0 = 15.0
  rmax = lda0**2 / dlda
  n = int(rmax/dr)

  q = np.linspace(0, 2*np.pi/dr, n)
  S = np.zeros((n, npairs+1))
  S[:, 0] = q

  for i in range(npairs):
    #Integrand in the fourier transform
    ker = (gr[:, 2*i + 1] - 1) * r
    #Imaginary (sin) part of the Fourier transform
    ft = np.imag(np.fft.fft(ker, n)) * dr
    #Structure factor
    #We split the q = 0 case, since it is ill-defined
    S[1:, i+1] = 1 - (ft[1:] / q[1:]) * (4 * np.pi * _d)

  mask = ((q > 0.2) & (q < 0.8))
  s_pasta = S[mask, 2]
  q_pasta = S[mask, 0]
  idx = np.argmax(s_pasta)
  s_max = s_pasta[idx]
  q_max = q_pasta[idx]

  return S, q_max, s_max

def fit(gr):
  """
  Fit the very-long correlations of neutron-neutron (near 15 fm)
  with a sine and returns the height of the fit and the wavelength
  """
  g = gr[:, 7]
  r = gr[:, 0]
  #Filter out initial zeros
  for i in range(len(g)):
    if g[i] > 0:
        break
    init = i
  g = g[init:]
  r = r[init:]
  sig = (g-1) * r
  sm_r, sm_sig = smooth(sig, r, dr=3.0)
  func = lambda x, A, l, ph: A * np.sin((2.0*np.pi*x) / l + ph)
  try:
    sm_popt, sm_pcov = curve_fit(func, sm_r, sm_sig, p0=[3.0, 15.0, 0.0])
    popt, pcov = curve_fit(func, r, sig, p0=sm_popt)
  except RuntimeError:
    return float("nan"), float("nan"), float("nan"), float("nan")
  tval = t.ppf(1.0-0.34/2, 3)
  h = sm_popt[0]
  dh = sm_pcov[0][0]**0.5 * tval
  l = sm_popt[1]
  dl = sm_pcov[1][1]**0.5 * tval
  return h, dh, l, dl


def mste(lmp, N):
  """
  Wrapper for the compute/mste method in LAMMPS (shipped with the package).
  LAMMPS returns simply an array with the cluster ID of each
  tag, so inside we do some calculations to return the mass
  distribution. The compute name *is* and *must be* mste [although
  this is meant to be obscured to the user in the setup() method].
  There is also a problem with setting the cutoff radius. It cannot
  be set from within the mste method, but it isn't very important: it
  should always be the cutoff of the pandha potential. Anyway, this is
  very unstable, so proceed with care when handling.

  There is a resize with mste, to make sure we can add different
  arrays afterwards. We add a lot of zeros and lose the sparsity
  of the mste, but early optimization...

  * system: system to analyze
  """

  ext = lmp.extract_compute("mste", 1, 1)
  tmp = np.fromiter(ext, dtype=np.int, count=N)
  ext = lmp.extract_atom("type", 0)
  typ = np.fromiter(ext, dtype=np.int, count=N)
  # We loop over all the different clusterIDs present, get the mass as
  # how often it occurs and the fraction as

  single_ids = set([i for i in tmp])
  clust = np.zeros((N,2), dtype=np.float)
  for (idx, clusterid) in enumerate(single_ids):
    indices = (tmp == clusterid)
    clus_type = typ[indices]
    mass = np.size(clus_type)
    frac = float(sum(typ[indices] - 1)) / mass
    clust[idx, 0] = mass
    clust[idx, 1] = frac

  clust = clust[clust[:, 0] != 0, :]
  mean = np.mean(clust)
  std = np.std(clust)
  single_masses = set([i for i in clust[:, 0]])
  histo_occ = np.zeros((N + 1), dtype=np.float)
  histo_frac = np.zeros((N + 1), dtype=np.float)
  for (idx, mass) in enumerate(single_masses):
    if mass == 0: continue
    indices = (clust[:, 0] == mass)

    occ = sum(indices)
    frac = sum(clust[indices, 1]) / occ
    histo_occ[mass] = occ
    histo_frac[mass] = frac

  return histo_occ, histo_frac, mean, std

def minkowski(lmp, rad, rcell):
  """
  Wrapper for minkowski. Need to add a warning with the
  size of the lattice for digitalization!
  """

  tmp = (ct.c_double * 4)()
  analysis.minkowski.argtypes = [ct.c_void_p, ct.c_double,
                                 ct.c_double, ct.c_void_p]
  analysis.minkowski(lmp.lmp, rad, rcell, tmp)
  return np.fromiter(tmp, dtype=np.float, count=4)

def lind(lmp):
    raise AttributeError("Don't know what to do with Lindemann :(")

def thermo(lmp, N):
  """
  Wrapper to LAMMPS internal computes.
  To avoid adding unnecesary computes to LAMMPS, we just reference
  to the default computes created for the LAMMPS inner thermo output.

  We have an advantage here: every time LAMMPS ends a run,
  calculates again thermo_temp, etc if they are in the thermo_style
  """

  temp = lmp.extract_compute("thermo_temp", 0, 0)
  epair = lmp.extract_compute("thermo_pe", 0, 0)/N
  press = lmp.extract_compute("thermo_press", 0, 0)
  ke = 3.0/2.0 * temp
  etot = epair + ke
  return temp, ke, epair, etot, press


def smooth(sig, r, dr=3.0):
  """
  Smooth signals
  """
  win_size = int(dr/(r[1] - r[0]))
  if win_size % 2 == 0:
    win_size += 1
  init = (win_size -1)/2
  fin_size = len(r) - win_size
  sm_sig = np.zeros(fin_size)
  sm_r = np.zeros(fin_size)
  for i in range(fin_size):
    for j in range(win_size):
      sm_sig[i] += sig[i+j-init]
    sm_sig[i] /= win_size
    sm_r[i] = r[i+init]
  return sm_r, sm_sig
