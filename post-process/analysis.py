"""
This files has analysis routines. They have to be in Python,
because we want to tamper with them
"""

import numpy as np
import itertools as it
import ctypes as C

libanalysis = C.CDLL('libanalysis.so')
rdf_c = libanalysis.rdf
ssf_c = libanalysis.ssf
cluster_c = libanalysis.cluster

def rdf(x, t, size):
  """
  rdf that gets a box and calculates with PBC in 3d.

  Parameters:

  - x: 2D numpy array
       First dimension is the number of particles, second one is x y and z

  - t: numpy array
       Type of particles
  """
  nbins = 200
  tmp = (C.c_double * (nbins * 5))()
  natoms = np.shape(x)[0]
  x_p = x.ctypes.data_as(C.POINTER(C.c_double))
  t_p = t.ctypes.data_as(C.POINTER(C.c_int))
  rdf_c.argtypes = [C.POINTER(C.c_double), C.POINTER(C.c_int),
                    C.c_int, C.c_int, C.c_double, C.POINTER(C.c_double)]
  rdf_c(x_p, t_p, natoms, nbins, size, tmp)
  gr = np.frombuffer(tmp, dtype = np.double, count = nbins * 5)
  return gr.reshape((nbins, 5))

def cluster(x, v, t, size):
  """
  cluster that gets a box and calculates MST or MSTE

  Parameters:

  - x: 2D numpy array
       First dimension is the number of particles, second one is x y and z

  - v: 2D numpy array
       First dimension is the number of particles, second is vx, vy and vz

  - t: numpy array
       Type of particles
  """
  _t = {'x': C.POINTER(C.c_double),
        'v': C.POINTER(C.c_double),
        'type': C.POINTER(C.c_int),
        'natoms': C.c_int,
        'size': C.c_double,
        'index': C.POINTER(C.c_int)}
  natoms = np.shape(x)[0]
  tmp = (C.c_int * (natoms))()
  x_p = x.ctypes.data_as(C.POINTER(C.c_double))
  v_p = v.ctypes.data_as(C.POINTER(C.c_double))
  t_p = t.ctypes.data_as(C.POINTER(C.c_int))
  cluster_c.argtypes = [_t['x'], _t['v'], _t['type'], _t['natoms'],
                        _t['size'], _t['index']]
  rdf_c(x_p, v_p, t_p, natoms, size, tmp)
  index = np.frombuffer(tmp, dtype=np.int, count=natoms)
  return index

def ssf(x, t, k, size, nrep=0):
  """
  ssf that gets a box and calculates with PBC in 3d.

  Parameters:

  - x: 2D numpy array
       First dimension is the number of particles, second one is x y and z

  - t: numpy array
       Type of particles

  - k: 1D numpy array
       Wavenumbers to calculate
  """
  npoints = np.shape(k)[0]
  out = (C.c_double * (npoints * 5))()
  natoms = np.shape(x)[0]
  x_p = x.ctypes.data_as(C.POINTER(C.c_double))
  t_p = t.ctypes.data_as(C.POINTER(C.c_int))
  k_p = k.ctypes.data_as(C.POINTER(C.c_double))
  ssf_c.argtypes = [C.POINTER(C.c_double), C.POINTER(C.c_int), C.c_int,
                    C.c_double, C.c_int, C.c_int, C.POINTER(C.c_double),
                    C.POINTER(C.c_double)]
  ssf_c(x_p, t_p, natoms, size, npoints, nrep, k_p, out)
  ssf = np.frombuffer(out, dtype=np.double, count=npoints * 5)
  return ssf.reshape((npoints, 5))


def transf_rdf(gr, density):
  """
  Calculate structure factor given the radial distribution function.
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
  ker = (gr[:, 1] - 1) * r
  #Imaginary (sin) part of the Fourier transform
  ft = np.imag(np.fft.fft(ker, n)) * dr
  #We split the q = 0 case, since it is ill-defined
  S[1:, 1] = 1 - (ft[1:] / q[1:]) * (4 * np.pi * _d)
  S[0, 1] = 1 + dr * sum(ker * r) * (4 * np.pi * _d)
  return S


def int_s(s, q, r):
  dq = q[1] - q[0]
  ker = np.sum(q*(s-1)*np.sin(q*r))*dq
  g = 2/np.pi*ker
  #g = ker
  g = g/(4*np.pi*r*0.05)+1
  return g


if __name__ == '__main__':
  import extract as E
  ext = E.Extraction('/home/pablo/artemis/eos/data')
  parameters = {'x': 0.5, 'd': 0.08, 'T': 0.5, 'V': 'horowitz'}
  x = ext.particle('x', parameters)
  t = ext.particle('type', parameters)
  v = ext.particle('v', parameters)
  size = (11000/0.08)**(1.0/3.0)
  d = np.shape(x)[0]/(size ** 3)
  gr = rdf(x, t, size)
  sk_transf = transf_rdf(gr, d)
  mst = cluster(x, v, t, size)
  k1 = np.linspace(0.2, 0.5, 101)
  k2 = np.linspace(3.0, 4.0, 101)
  k = np.hstack((k1, k2))

  sk = ssf(x, t, k, size, 1)
