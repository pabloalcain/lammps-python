"""
RDF calculation
"""

import ctypes as ct
import numpy as np
import os

_DIRNAME = os.path.dirname(__file__)
librdf = ct.CDLL(os.path.join(_DIRNAME, 'librdf.so'))
rdf_c = librdf.rdf

def rdf(x, t, size, pairs, nbins, pbc=True):
  """Calculate RDF.

  Parameters
  ----------

  x : numpy array
      Positions of the particles in the system

  t : numpy array
      Types of the particles

  size : float
      Length of the box

  pairs: iterable
      List of pairs to consider. Each element of the list is a tuple
      of size 2, with the first element being the types considered
      and the second one the other types. For example, if we want to
      calculate the RDF on a system with 2 types of particles of all
      types vs all types and type 1 vs type 2, the pairs should be:

      pairs = [((1, 2), (1, 2)), ((1,), (2,))]

      A keyword for "all types" is type 0

  nbins : int
      Number of bins to consider

  pbc : boolean, optional
      Whether or not to use periodic boundary conditions. Default is
      True.

  Returns
  -------

  rdf : numpy array
      An array with the information of the compute. The first column
      is 'r' and the rest are the RDF calculated for the pair list.
  """
  #TODO: Add some warning in case the indices repeat
  natoms = np.shape(x)[0]
  npairs = len(pairs)
  ncol = npairs + 1
  pair_ar = np.zeros(2*npairs)
  i = 0
  for p in pairs:
    pair_ar[i] = p[0][0]
    pair_ar[i+1] = p[1][0]
    i += 2
  tmp = (ct.c_double * (nbins * ncol))()
  x_p = x.ctypes.data_as(ct.c_void_p)
  t_p = t.ctypes.data_as(ct.c_void_p)
  pair_p = pair_ar.ctypes.data_as(ct.c_void_p)
  rdf_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int, ct.c_int,
                    ct.c_double, ct.c_void_p, ct.c_int, ct.c_bool,
                    ct.c_void_p]
  rdf_c(x_p, t_p, natoms, nbins, size, pair_p, npairs, pbc, tmp)
  value = np.frombuffer(tmp, dtype=np.double, count=nbins * ncol)
  return value.reshape((nbins, ncol))


def ssf(gr, density, pbc=True):
  """
  Calculate structure factor given the radial distribution function.

  Parameters
  ----------

  gr : numpy array
      Radial Distribution Function. First column is the position,
      second is the RDF itself

  density : float
      Density of the studied system

  pbc : boolean, optional
      Whether or not the rdf was calculated with periodic boundary
      conditions.

  Returns
  -------

  S : numpy array
      The structure factor of the system. First column is the
      wavenumber, second is the ssf itself.
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
  S = np.zeros((n, 2))
  S[:, 0] = q

  #Integrand in the fourier transform
  if pbc: mean = 1
  else: mean = 0
  ker = (gr[:, 1] - mean) * r
  #Imaginary (sin) part of the Fourier transform
  ft = np.imag(np.fft.fft(ker, n)) * dr
  #We split the q = 0 case, since it is ill-defined
  S[1:, 1] = 1 - (ft[1:] / q[1:]) * (4 * np.pi * _d)
  S[0, 1] = 1 + dr * sum(ker * r) * (4 * np.pi * _d)
  return S
