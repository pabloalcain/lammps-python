"""
MST Compute class
"""

from Compute import Compute
import ctypes as ct
import numpy as np


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


class MST(Compute):
  """
  MST/MSTE calculation.
  """
  def __init__(self, energy=True, pbc=True):
    """
    Constructor.

    Parameters
    ----------

    energy : boolean
        Whether or not to use energy considerations (MSTE) in the
        cluster recognition. Default is True.
    """
    self.energy = energy
    self.pbc = pbc
    super(MST, self).__init__()

  def compute(self, system):
    """Calculate MST.

    Parameters
    ----------

    system : System
        System on which we calculate the Minimum Spanning Tree

    Returns
    -------

    value, (mst, inf) : numpy array, numpy array, list
        value is the [mass, fraction] reduction
        mst is the array of indices to which each particle belongs.
        inf is the list of infinite clusters.

    Notes
    -----

    So far, this only works for systems with `medium` potential.

    """
    natoms = system['N']
    expansion = system['expansion']
    size = system.size
    index_p = (ct.c_int * natoms)()
    x_p = system.x.ctypes.data_as(ct.c_void_p)
    v_p = system.v.ctypes.data_as(ct.c_void_p)
    t_p = system.t.ctypes.data_as(ct.c_void_p)
    cluster_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p,
                          ct.c_int, ct.c_double, ct.c_bool,
                          ct.c_void_p]
    cluster_c(x_p, v_p, t_p, natoms, size, self.energy, index_p)
    index = np.frombuffer(index_p, dtype=ct.c_int, count=natoms)
    nclus = len(np.unique(index))
    #TODO: check for size not too large
    if nclus > 1000: nclus = 1000
    guess = nclus ** 2 * 6
    tmp = (ct.c_int * (3 * guess))()

    connections_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p,
                              ct.c_void_p, ct.c_int, ct.c_float,
                              ct.c_bool, ct.c_double, ct.c_void_p]
    connections_c.restype = ct.c_int
    count = connections_c(index_p, x_p, v_p, t_p, natoms, size,
                          expansion, self.energy, tmp)
    conn = np.frombuffer(tmp, dtype=ct.c_int, count=count * 3)
    conn = conn.reshape((count, 3))
    graph, cnct = _create_graph(conn)
    mst = index.copy()
    idx = 0
    inf = []
    value = []
    for nodes in _partition(graph):
      idx += 1
      this_graph = {}
      for n in nodes:
        this_graph[n] = graph[n]
      _, inf_clusters = _find_paths(this_graph, cnct)
      mass = 0
      protons = 0
      for n in nodes:
        mst[index == n] = idx
        mass += len(index == n)
        protons += len((index == n) & (system.t == 2))
      frac = protons/mass
      if len(inf_clusters) != 0:
        inf.append(idx)
        mass = np.inf
      value.append([mass, frac])
    return value, (mst, inf)

  def tally(self, value):
    """
    We need to override the previous tally, since this compute does
    not average trivially. We create the histogram from the values.
    """
    self.idx += 1
    for item in self.value:
      item[0] *= (self.idx - 1)/self.idx
    for mass, fraction in value:
      self.value[mass, 0] += 1/self.idx
      self.value[mass, 1] *= (self.value[mass] - 1)/self.value[mass]
      self.value[mass, 1] += fraction/self.value[mass]
