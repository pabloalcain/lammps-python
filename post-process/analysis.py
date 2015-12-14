"""
This files has analysis routines. They have to be in Python, because we want to tamper with them
"""

import numpy as np
import itertools as it
import ctypes as C

libanalysis = C.CDLL('libanalysis.so')
rdf_c = libanalysis.rdf
ssf_c = libanalysis.ssf

def _f1(l):
  return np.arctan(np.sqrt(4*l**2 - 2))

def _f2(l):
  ang = 2*l*(4*l**2-3)/(np.sqrt(4*l**2-2)*(4*l**2+1))
  return 8*l*np.arctan(ang)

def _vol(l):
  if l < 0.5:
    v = 4.0/3.0 * np.pi * l**3
  elif l < 1.0/np.sqrt(2):
    v = -np.pi/12*(3-36*l**2+32*l**3)
  elif l < np.sqrt(3)/2:
    v = -np.pi/4 + 3 * np.pi * l**2
    v += np.sqrt(4 * l**2 - 2)
    v += (1 - 12 * l**2) * _f1(l)
    v += 2.0/3.0 * l**2 * _f2(l)
  return v

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

def ssf(x, t, k, size):
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
                    C.c_double, C.c_int, C.POINTER(C.c_double),
                    C.POINTER(C.c_double)]
  ssf_c(x_p, t_p, natoms, size, npoints, k_p, out)
  ssf = np.frombuffer(out, dtype=np.double, count=npoints * 5)
  return ssf.reshape((npoints, 5))


def rdf_py(x, box):
  """
  rdf that gets a box and calculates with PBC in 3d.

  Parameters:

  - x: 2D numpy array
  First dimension is the number of particles, second one is x y and z
  """

  if type(box) == float:
    xlo = -box/2
    ylo = -box/2
    zlo = -box/2
    xhi = box/2
    yhi = box/2
    zhi = box/2
  else:
    raise ValueError("Only support cubic boxes")

  vol = (xhi - xlo) * (yhi - ylo) * (zhi - zlo)
  rmax = np.sqrt(3)*box/2
  npoints = 200
  delr = rmax/(npoints+1)
  rdf = np.zeros((npoints, 2))
  rdf[:, 0] = np.linspace(0, rmax, npoints + 1)[:-1]
  natoms = np.shape(x)[0]
  for i in range(natoms):
    for j in range(i+1, natoms):
      dx = x[j] - x[i]
      if dx[0] > box/2: dx[0] -= box
      elif dx[0] < -box/2: dx[0] += box
      if dx[1] > box/2: dx[1] -= box
      elif dx[1] < -box/2: dx[1] += box
      if dx[2] > box/2: dx[2] -= box
      elif dx[2] < -box/2: dx[2] += box
  r = np.linalg.norm(dx)
  idx = int(r*npoints/rmax)
  rdf[idx, 1] += 1
  density = natoms / vol
  for i in rdf:
    l = i[0]/box
    p = _vol(l + delr/box) - _vol(l)
    i[1] = i[1] / (p * density * box**3 * (natoms/2))
  return rdf


def pbc_rdf(x, box):
  """
  rdf that gets a box and calculates with PBC in 3d, replicating as needed

  Parameters:

  - x: 2D numpy array
  First dimension is the number of particles, second one is x y and z
  """
  if type(box) == float:
    xlo = -box/2
    ylo = -box/2
    zlo = -box/2
    xhi = box/2
    yhi = box/2
    zhi = box/2
  else:
    raise ValueError("Only support cubic boxes")

  natoms = np.shape(x)[0]
  base = (box * np.linspace(0, 2, 3))**2
  rmax = 0.5*box
  vol = rmax * rmax * rmax
  npoints = 1000
  delr = rmax/npoints
  rdf = np.zeros((npoints, 2))
  rdf[:, 0] = np.linspace(0, rmax, npoints)
  hits = 0
  for i in range(natoms):
    for j in range(natoms):
      for ix, iy, iz in it.product(range(-1, 2), range(-1, 2), range(-1, 2)):
        if i == j and ix == 0 and iy == 0 and iz == 0: continue
        d2 = base[abs(ix)] + base[abs(iy)] + base[abs(iz)]
        dx = x[j, 0] - x[i, 0]
        dy = x[j, 1] - x[i, 1]
        dz = x[j, 2] - x[i, 2]
        dd = ix * dx + iy * dy + iz * dz
        dd *= 2*box
        dr = x[j] - x[i]
        r = np.sqrt(np.dot(dr, dr) + d2 + dd)
        idx = int(r*npoints/rmax)
        if r < rmax:
          hits += 1
          rdf[idx, 1] += 1

  c = 4*np.pi / (3 * vol)
  for i in rdf:
    dV = (i[0] + delr) ** 3 - i[0]**3
    i[1] = i[1]  / (dV * c * hits)
  return rdf

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
  #n = len(r)
  q = np.linspace(0, 2*np.pi/dr, n)
  S = np.zeros((n, 2))
  S[:, 0] = q

  #Integrand in the fourier transform
  ker = (gr[:, 1] - 1) * r
  #Imaginary (sin) part of the Fourier transform
  ft = np.imag(np.fft.fft(ker, n)) * dr
  #We split the q = 0 case, since it is ill-defined
  S[1:, 1] = 1 - (ft[1:] / q[1:]) * (4 * np.pi * _d)
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
  x = E.positions('dump.lammpstrj', 0)
  t = E.types('dump.lammpstrj', 0)
  gr = rdf(x, t, 23.9397*2)
  print "Done rdf!"
  k1 = np.linspace(0.0, 0.5, 100)
  k2 = np.linspace(2.0, 4.0, 100)
  k = k1 #np.hstack((k1, k2))
  #k = np.linspace(0.0, 0.5, 100)
  d = 0.05
  N = 5488
  l = (N/d)**(1.0/3.0)
  sk_transf = transf_rdf(gr, d)
  k =sk_transf[:, 0].copy()
  sk = ssf(x, t, k, l)
