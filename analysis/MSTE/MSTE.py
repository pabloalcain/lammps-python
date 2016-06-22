"""
MSTE calculation
"""

import ctypes as ct
import numpy as np
import os
import warnings

_DIRNAME = os.path.dirname(__file__)
libmste = ct.CDLL(os.path.join(_DIRNAME, 'libmste.so'))
connections_c = libmste.connections
cluster_c = libmste.cluster


#TODO: Maybe avoid this int2str and str2int conversion?
def _idx2wall(wall):
  """
  Give a value in binary +/- 1/10/100 to each face
  """
  if wall == 1:
    return 'W'
  if wall == -1:
    return 'E'
  if wall == 2:
    return 'N'
  if wall == -2:
    return 'S'
  if wall == 4:
    return 'U'
  if wall == -4:
    return 'D'
  else:
    raise ValueError('Index {0} not understood'.format(wall))

def _wall2idx(wall):
  """
  Give a value in binary +/- 1/10/100 to each face
  """
  if wall == 'W':
    return 1
  if wall == 'E':
    return -1
  if wall == 'N':
    return 2
  if wall == 'S':
    return -2
  if wall == 'U':
    return 4
  if wall == 'D':
    return -4
  else:
    raise ValueError('Wall {0} not understood'.format(wall))

def _create_graph(conn):
  """
  From a connection list, create a graph dictionary and format
  connections as a dict

  Parameters
  ----------

  conn : 2D numpy array
      A connection array with (i, j, wall) indices in each entry

  Returns
  -------

  graph, cnct: dict, dict
      Graph and connections in a dictionary fashion
  """
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

def _plain_bfs(graph, source):
  """A fast BFS node generator"""
  seen = set()
  nextlevel = {source}
  while nextlevel:
    thislevel = nextlevel
    nextlevel = set()
    for vertex in thislevel:
      if vertex not in seen:
        yield vertex
        seen.add(vertex)
        nextlevel.update(graph[vertex])

def _partition(graph):
  """
  Partition a graph in disconnected subgraphs
  """
  seen = set()
  for vertex in graph:
    if vertex not in seen:
      conn = set(_plain_bfs(graph, vertex))
      yield conn
      seen.update(conn)

def _find_cycles(graph):
  """
  Get cycles in undirected graph. But we return also every edge as a
  "cycle", since it can generate an infinite cluster
  """
  gnodes = set(graph.keys())
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
      for nbr in graph[z]:
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
  for gnode in graph.keys():
    for i in graph[gnode]:
      if i > gnode:
        cycs.append([gnode, i])
  return cycs

def _find_paths(graph, cncts):
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
  for cyc in _find_cycles(graph):
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

    if inf:
      inf_clusters.append(cyc)
  return cycles, inf_clusters

def mste(x, v, t, box, energy, expansion=0.0):
  """Calculate MSTE.

  Parameters
  ----------

  x : numpy float64 array
      Positions of the particles in the system

  v : numpy float64 array
      Velocities of the particles

  t : numpy int32 array
      Types of the particles

  box : numpy float64 array
      Box

  energy : boolean
      Whether to do energy considerations in the cluster computation.

  expansion : float, optional
      The expansion velocity of the walls of the box, in box units.

  Returns
  -------

  value, (mst, inf) : numpy array, numpy array, list
      value is the [mass, occupancy, fraction] histogram
      mst is the array of indices to which each particle belongs.
      inf is the list of infinite clusters.

  Notes
  -----

  1. So far, when energy is "True", this only works for systems with
  `medium` potential.

  2. The box needs to be cubic for this to work
  """
  if x.dtype != np.float64:
    warnings.warn("Type of x should be 64bit float", UserWarning)
    x = x.astype(np.float64)
  if v.dtype != np.float64:
    warnings.warn("Type of v should be 64bit float", UserWarning)
    v = v.astype(np.float64)
  if t.dtype != np.int32:
    warnings.warn("Type of t should be 32bit float", UserWarning)
    t1 = t.astype(np.int32)
    for i, j in zip(t.flatten(), t1.flatten()):
      if i != j:
        raise ValueError("Casting resulted in different values")

  index = cluster(x, v, t, box, energy)
  conn = connections(x, v, t, index, box, energy, expansion)
  graph, cnct = _create_graph(conn)
  #TODO: This should be an independent function
  mst = index.copy()
  natoms = np.shape(x)[0]
  inf = []
  value = np.zeros((natoms + 1, 3))
  value[:, 0] = range(natoms + 1)
  for nodes in _partition(graph):
    this_graph = {}
    for n in nodes:
      this_graph[n] = graph[n]
      mst[index == n] = min(nodes)
    _, inf_clusters = _find_paths(this_graph, cnct)
    if len(inf_clusters) != 0:
      inf.append(min(nodes))
  i = 0
  for clus in np.unique(mst):
    mass = np.count_nonzero(mst == clus)
    protons = np.count_nonzero((mst == clus) & (t.flatten() == 2))
    frac = float(protons)/mass
    if clus in inf:
      mass = 0
    value[mass, 1] += 1
    value[mass, 2] *= (value[mass, 1] - 1.0)/value[mass, 1]
    value[mass, 2] += frac/value[mass, 1]

  #Reduce to sequential cluster numbering
  mst2 = mst.copy()
  i = 0
  for v in np.unique(mst2):
    mst[mst2==v] = i
    i += 1
  return value, (mst, inf)

def cluster(x, v, t, box, energy):
  """
  Get either MST or MSTE clusters from the configuration.

  Parameters
  ----------

  x : numpy float64 array
      Positions of the particles in the system

  v : numpy float64 array
      Velocities of the particles

  t : numpy int32 array
      Types of the particles

  box : numpy float64 array
      Box

  energy : boolean
      Whether to do energy considerations in the cluster computation.

  Returns
  -------

  index : numpy array
      An array with the cluster index of each particle.
  """
  size_x = box[0][1] - box[0][0]
  size_y = box[1][1] - box[1][0]
  size_z = box[2][1] - box[2][0]
  if size_x != size_y or size_y != size_z:
    raise ValueError("The box should be cubic for this to work")
  else:
    size = size_x

  natoms = np.shape(x)[0]
  index_p = (ct.c_int * natoms)()
  x_p = x.ctypes.data_as(ct.c_void_p)
  v_p = v.ctypes.data_as(ct.c_void_p)
  t_p = t.ctypes.data_as(ct.c_void_p)
  cluster_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p,
                        ct.c_int, ct.c_double, ct.c_bool,
                        ct.c_void_p]
  cluster_c(x_p, v_p, t_p, natoms, size, energy, index_p)
  index = np.frombuffer(index_p, dtype=ct.c_int, count=natoms)
  return index


def connections(x, v, t, index, box, energy, expansion):
  """Find the connections with a given cluster distribution.

  Parameters
  ----------

  x : numpy float64 array
      Positions of the particles in the system

  v : numpy float64 array
      Velocities of the particles

  t : numpy int32 array
      Types of the particles

  index : numpy array
      Cluster index of the particles

  box : numpy float64 array
      Box

  energy : boolean
      Whether to do energy considerations in the cluster computation.

  expansion : float, optional
      The expansion velocity of the walls of the box, in box units.

  Returns
  -------

  conn : 2D numpy array
      The connections between the clusters, in the format:
      [[cluster1, cluster2, wall], [cluster1, cluster2, wall]...]

  Notes
  -----

  1. So far, when energy is "True", this only works for systems with
     `medium` potential.

  2. The box needs to be cubic for this to work
  """
  size_x = box[0][1] - box[0][0]
  size_y = box[1][1] - box[1][0]
  size_z = box[2][1] - box[2][0]
  if size_x != size_y or size_y != size_z:
    raise ValueError("The box should be cubic for this to work")
  else:
    size = size_x

  natoms = np.shape(x)[0]
  index_p = index.ctypes.data_as(ct.c_void_p)
  x_p = x.ctypes.data_as(ct.c_void_p)
  v_p = v.ctypes.data_as(ct.c_void_p)
  t_p = t.ctypes.data_as(ct.c_void_p)
  nclus = len(np.unique(index))
  #TODO: check for size not too large
  if nclus > 1000:
    nclus = 1000
  guess = nclus ** 2 * 8
  tmp = (ct.c_int * (3 * guess))()
  connections_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p,
                            ct.c_void_p, ct.c_int, ct.c_double,
                            ct.c_double, ct.c_bool, ct.c_void_p]
  connections_c.restype = ct.c_int



  count = connections_c(index_p, x_p, v_p, t_p, natoms, size,
                        expansion, energy, tmp)
  conn = np.frombuffer(tmp, dtype=ct.c_int, count=count * 3)
  conn = conn.reshape((count, 3))
  return conn
