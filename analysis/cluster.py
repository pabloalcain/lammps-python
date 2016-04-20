"""
Routines for cluster analysis and graph reduction
"""

import numpy as np
from analysis import cluster_c, connections_c
import ctypes as ct

def _idx2wall(w):
  """
  Give a value in binary +/- 1/10/100 to each face
  """
  if w == 1: return 'W'
  elif w == -1: return 'E'
  elif w == 2: return 'N'
  elif w == -2: return 'S'
  elif w == 4: return 'U'
  elif w == -4: return 'D'

def _wall2idx(w):
  """
  Give a value in binary +/- 1/10/100 to each face
  """
  if w == 'W': return 1
  elif w == 'E': return -1
  elif w == 'N': return 2
  elif w == 'S': return -2
  elif w == 'U': return 4
  elif w == 'D': return -4


def reduce_mste(index, t):
  """
  Get a list of indices and types and returns the histogram of
  occupation and of fractions
  """
  N = np.shape(index)[0]
  single_ids = set([i for i in index])
  clust = np.zeros((N, 2), dtype=np.float)

  for (idx, clusterid) in enumerate(single_ids):
    indices = (index == clusterid)
    clus_type = t[indices]
    mass = np.size(clus_type)
    frac = float(sum(t[indices] - 1)) / mass
    clust[idx, 0] = mass
    clust[idx, 1] = frac
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
  indices = (histo_occ != 0.0)
  mass = np.arange(len(histo_occ))
  frac = histo_frac[indices].copy()
  mass = mass[indices].copy()
  occ = histo_occ[indices].copy() # Filter out
  return np.vstack((mass, occ, frac)).T

def cluster(x, v, t, size, energy=False, pbc=False):
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
  _t = {'x': ct.POINTER(ct.c_double),
        'v': ct.POINTER(ct.c_double),
        'type': ct.POINTER(ct.c_int),
        'natoms': ct.c_int,
        'size': ct.c_double,
        'index': ct.POINTER(ct.c_int),
        'energy': ct.c_bool,
        'pbc': ct.c_bool}
  natoms = np.shape(x)[0]
  tmp = (ct.c_int * (natoms))()
  x_p = x.ctypes.data_as(ct.POINTER(ct.c_double))
  v_p = v.ctypes.data_as(ct.POINTER(ct.c_double))
  t_p = t.ctypes.data_as(ct.POINTER(ct.c_int))
  cluster_c.argtypes = [_t['x'], _t['v'], _t['type'], _t['natoms'],
                        _t['size'], _t['energy'], _t['pbc'], _t['index']]
  cluster_c(x_p, v_p, t_p, natoms, size, energy, pbc, tmp)
  index = np.frombuffer(tmp, dtype=ct.c_int, count=natoms)
  return index

def connections(index, x, v, t, size, expansion=0.0, energy=False):
  """
  Get a list of clusters and a box and finds the connection through
  the walls.

  Parameters:

  - index: numpy array
       Cluster index of the particle

  - x: 2D numpy array
       First dimension is the number of particles, second one is x y and z

  - v: 2D numpy array
       First dimension is the number of particles, second is vx, vy and vz

  - t: numpy array
       Type of particles
  """

  _t = {'x': ct.POINTER(ct.c_double),
        'v': ct.POINTER(ct.c_double),
        'type': ct.POINTER(ct.c_int),
        'natoms': ct.c_int,
        'size': ct.c_double,
        'index': ct.POINTER(ct.c_int),
        'expansion': ct.c_double,
        'energy': ct.c_bool,
        'connect': ct.POINTER(ct.c_int)}

  natoms = np.shape(x)[0]
  nclus = len(np.unique(index))
  #TODO: check for size not too large
  if nclus > 1000: nclus = 1000
  guess = nclus ** 2 * 6
  tmp = (ct.c_int * (3 * guess))()
  x_p = x.ctypes.data_as(ct.POINTER(ct.c_double))
  v_p = v.ctypes.data_as(ct.POINTER(ct.c_double))
  t_p = t.ctypes.data_as(ct.POINTER(ct.c_int))
  index_p = index.ctypes.data_as(ct.POINTER(ct.c_int))
  connections_c.argtypes = [_t['index'], _t['x'], _t['v'], _t['type'],
                            _t['natoms'], _t['size'], _t['expansion'],
                            _t['energy'], _t['connect']]
  connections_c.restype = ct.c_int
  count = connections_c(index_p, x_p, v_p, t_p,
                        natoms, size, expansion, energy, tmp)
  conn = np.frombuffer(tmp, dtype=ct.c_int, count=count * 3)
  conn = conn.reshape((count, 3))
  graph = {}
  cnct = {}
  for i in conn:
    iclus = i[0]
    jclus = i[1]
    wall = _idx2wall(i[2])
    try:
      if not jclus in graph[iclus]:
        graph[iclus].append(jclus)
    except KeyError:
      graph[iclus] = [jclus]
    try:
      if not wall in cnct[(iclus, jclus)]:
        cnct[(iclus, jclus)].append(wall)
    except KeyError:
      cnct[(iclus, jclus)] = [wall]
  return graph, cnct


def _plain_bfs(G, source):
  """A fast BFS node generator"""
  seen = set()
  nextlevel = {source}
  while nextlevel:
    thislevel = nextlevel
    nextlevel = set()
    for v in thislevel:
      if v not in seen:
        yield v
        seen.add(v)
        nextlevel.update(G[v])


def partition(g):
  """
  Partition a graph in disconnected subgraphs
  """
  seen = set()
  for v in g:
    if v not in seen:
      c = set(_plain_bfs(g, v))
      yield c
      seen.update(c)

def find_cycles(g):
  """
  Get cycles in undirected graph. But we return also every edge as a
  "cycle", since it can generate an infinite cluster
  """
  gnodes = set(g.keys())
  cycs = []
  root = None
  while gnodes:  # loop over connected components
    if root is None:
      root = gnodes.pop()
    stack = [root]
    pred = {root:root}
    used = {root:set()}
    while stack:  # walk the spanning tree finding cycles
      z = stack.pop()  # use last-in so cycles easier to find
      zused = used[z]
      for nbr in g[z]:
        if nbr not in used:   # new node
          pred[nbr] = z
          stack.append(nbr)
          used[nbr] = set([z])
        elif nbr == z:    # self loops
          cycs.append([z])
        elif nbr not in zused:# found a cycle
          pn = used[nbr]
          cycle = [nbr, z]
          p = pred[z]
          while p not in pn:
            cycle.append(p)
            p = pred[p]
          cycle.append(p)
          cycs.append(cycle)
          used[nbr].add(z)
    gnodes -= set(pred)
    root = None
  for gnode in g.keys():
    for i in g[gnode]:
      if i > gnode:
        cycs.append([gnode, i])

  return cycs

def find_paths(g, cncts):
  """Find paths in an oriented [not directed!] graph and return the
  cycles and the paths through. Also inform whether they are infinite.
  """
  from copy import deepcopy
  def get_paths(cyc):
    """ Recursive function to get all the paths
    """
    paths = []
    def add_node(i=0, path=''):
      """
      Add new node
      """
      if i == len(cyc) - 1:
        paths.append(path)
        return
      for _ in cncts[(cyc[i], cyc[i+1])]:
        add_node(i + 1, deepcopy(path + _))
    add_node()
    return paths

  cycles = []
  inf_clusters = []
  for cyc in find_cycles(g):
    cyc.append(cyc[0])
    paths = get_paths(cyc)
    inf = False
    for path in paths:
      topo = 0
      for wall in ['N', 'S', 'E', 'W', 'U', 'D']:
        if wall in path:
          topo += _wall2idx(wall)
      inf = inf or topo != 0
      cycles.append((cyc, path))

    if inf: inf_clusters.append(cyc)
  return cycles, inf_clusters
