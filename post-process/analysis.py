"""
This files has analysis routines. They have to be in Python, because
we want to tamper with them
"""

import numpy as np
import itertools as it
import ctypes as C
import cluster

libanalysis = C.CDLL('libanalysis.so')
rdf_c = libanalysis.rdf
ssf_c = libanalysis.ssf

def rdf(x, t, size, nbins=200, pbc=True):
  """
  rdf that gets a box and calculates with PBC in 3d.

  Parameters:

  - x: 2D numpy array
       First dimension is the number of particles, second one is x y and z

  - t: numpy array
       Type of particles
  """
  tmp = (C.c_double * (nbins * 5))()
  natoms = np.shape(x)[0]
  x_p = x.ctypes.data_as(C.POINTER(C.c_double))
  t_p = t.ctypes.data_as(C.POINTER(C.c_int))
  rdf_c.argtypes = [C.POINTER(C.c_double), C.POINTER(C.c_int),
                    C.c_int, C.c_int, C.c_double, C.c_bool,
                    C.POINTER(C.c_double)]
  rdf_c(x_p, t_p, natoms, nbins, size, pbc, tmp)
  gr = np.frombuffer(tmp, dtype=np.double, count=nbins * 5)
  return gr.reshape((nbins, 5))


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


def transf_rdf(gr, density, pbc=True):
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
  if pbc: mean = 1
  else: mean = 0
  ker = (gr[:, 1] - mean) * r
  #Imaginary (sin) part of the Fourier transform
  ft = np.imag(np.fft.fft(ker, n)) * dr
  #We split the q = 0 case, since it is ill-defined
  S[1:, 1] = 1 - (ft[1:] / q[1:]) * (4 * np.pi * _d)
  S[0, 1] = 1 + dr * sum(ker * r) * (4 * np.pi * _d)
  mask = ((q > 0.2) & (q < 0.8))
  s_pasta = S[mask, 1]
  q_pasta = S[mask, 0]
  idx = np.argmax(s_pasta)
  s_max = s_pasta[idx]
  q_max = q_pasta[idx]
  return S, q_max, s_max

def int_s(s, q, r):
  dq = q[1] - q[0]
  ker = np.sum(q*(s-1)*np.sin(q*r))*dq
  g = 2/np.pi*ker
  #g = ker
  g = g/(4*np.pi*r*0.05)+1
  return g

if __name__ == '__main__':
  import extract as E
  import pylab as pl
  ext = E.Extraction('/home/palcain/opacity/data')
  parameters = {'x': 0.5, 'd': 0.02, 'T': 0.5, 'V': 'medium', 'N': 5488}
  for x in (0.2, 0.3, 0.4, 0.5):
    fs = open('s_{0}.dat'.format(x), 'w')
    fl = open('l_{0}.dat'.format(x), 'w')
    print>>fs, '#density {0}'.format(' '.join(map(str, np.linspace(0.5, 2.0, 31))))
    print>>fl, '#density {0}'.format(' '.join(map(str, np.linspace(0.5, 2.0, 31))))
    print x
    parameters['x'] = x
    for d in (0.001, 0.005, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08):
      print d
      parameters['d'] = d
      print>>fs, d,
      print>>fl, d,
      size = (5488/d)**(1.0/3.0)
      for T in np.linspace(0.5, 2.0, 31):
        print T
        parameters['T'] = T
        d = parameters['d']
        sk = 0
        q_max = 0
        s_max = 0
        for i in range(50):
          x = ext.particle('x', parameters, idx=i)
          t = ext.particle('type', parameters, idx=i)
          v = ext.particle('v', parameters, idx=i)
          gr_one = rdf(x, t, size, nbins=1000, pbc=False)
          sk_t, q_max_t, s_max_t = transf_rdf(gr_one[:, [0, 2]], d, pbc=False)
          sk += sk_t
          q_max += q_max_t
          s_max += s_max_t
        sk /= 50
        q_max /=50
        s_max /=50
        print>>fs, s_max,
        print>>fl, 2*np.pi/q_max,
      print>>fs
      print>>fl
    fs.close()
    fl.close()
  pl.figure()
  pl.plot(sk[:, 0], sk[:, 1], label='Without PBC')
  pl.vlines(q_max, 0, s_max*1.5)
  pl.legend()
 
  """
  mst = cluster.cluster(x, v, t, size, energy=True, pbc=False)
  graph, connections  = cluster.connections(mst, x, v, t, size, energy=True)
  #print graph
  for nodes in cluster.partition(graph):
    this_graph = {}
    for n in nodes:
      this_graph[n] = graph[n]
    cycles, inf_clusters = cluster.find_paths(this_graph, connections)
    mass = 0
    for n in nodes:
      mass += np.sum(mst == n)
    if len(inf_clusters) != 0:
      print "subgraph {0} is infinite of mass {1}".format(nodes, mass)
    else:
      print "subgraph {0} is finite of mass {1}".format(nodes, mass)

  #k1 = np.linspace(0.2, 0.5, 101)
  #k2 = np.linspace(3.0, 4.0, 101)
  #k = np.hstack((k1, k2))

  #sk = ssf(x, t, k, size, 1)
  """
