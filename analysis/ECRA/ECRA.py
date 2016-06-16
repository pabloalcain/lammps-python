"""
Cluster partition that minimizes the energy of the partition. Still
only for testing purposes. Still need to put C code into it to improve
its speed.
"""

import ctypes as ct
import numpy as np
import os
import warnings
import copy
from analysis import mste

_DIRNAME = os.path.dirname(__file__)
libecra = ct.CDLL(os.path.join(_DIRNAME, 'libecra.so'))
enpart_c = libecra.enpart

import random

def potential(dx, t1, t2):
  Vc = 1.44
  uc = 1.0/20.0
  (Vr, Va, V0) = (3088.118, 2666.647, 373.118)
  (ur, ua, u0) = (1.7468, 1.6, 1.5)
  if t1 == t2 == 1:
    e = V0 * np.exp(-u0 * dx) / dx
  elif t1 == t2 == 2:
    e = V0 * np.exp(-u0 * dx) / dx + Vc * np.exp(-uc * dx) / dx
  else:
    e = Vr * np.exp(-ur * dx) / dx - Va * np.exp(-ua * dx) / dx
  return e

def partition(collection):
  """
  An elegant iterative solution to the partitioning problem, from
  http://stackoverflow.com/questions/19368375/set-partitions-in-python
  """
  if len(collection) == 1:
    yield [collection]
    return

  first = collection[0]
  for smaller in partition(collection[1:]):
    # insert `first` in each of the subpartition's subsets
    for n, subset in enumerate(smaller):
      yield smaller[:n] + [[first] + subset]  + smaller[n+1:]
    # put `first` in its own subset
    yield [[first]] + smaller

def ecra_mst(x, v, t, box, expansion, index=None):
  """
  Calculate cluster recognition algorithm with Simulated Annealing as
  in [Dorso]_. We first calculate MST clusters and for those that are
  not infinite, we perform an ECRA calculation.

  .. [Dorso] Dorso and Randrup, Phys. Lett. B, 301, 4, 328-333

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

  expansion : float, optional
      The expansion velocity of the walls of the box, in box units.

  index : numpy int32 array, optional
      The original guess for the minimum energy

  Returns
  -------

  value, (ec, inf) : numpy array, numpy array, list
      value is the [mass, occupancy, fraction] histogram
      mst is the array of indices to which each particle belongs.
      inf is the list of infinite clusters, which is always empty.
  """
  value, (idx_mst, inf) = mste(x, v, t, box, False)
  _, (idx_mste, _) = mste(x, v, t, box, False)
  for clus in np.unique(idx_mst):
    if clus in inf: continue
    #For each different and finite MST cluster, we calculate ECRA
    mask = idx_mst == clus
    x_clus = x[mask].copy()
    v_clus = v[mask].copy()
    t_clus = t[mask].copy()
    i_clus = idx_mste[mask].copy()
    _, (idx, _) = ecra(x_clus, v_clus, t_clus, box, expansion,
                       index=i_clus)



  return value, (ec, [])


def ecra(x, v, t, box, expansion, index=None):
  """
  Calculate cluster recognition algorithm with Simulated Annealing as
  in [Dorso]_.

  .. [Dorso] Dorso and Randrup, Phys. Lett. B, 301, 4, 328-333

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

  expansion : float, optional
      The expansion velocity of the walls of the box, in box units.

  index : numpy int32 array, optional
      The original guess for the minimum energy

  Returns
  -------

  value, (ec, inf) : numpy array, numpy array, list
      value is the [mass, occupancy, fraction] histogram
      mst is the array of indices to which each particle belongs.
      inf is the list of infinite clusters, which is always empty.
  """
  npart = np.shape(x)[0]
  if index == None:
    index = np.arange(npart)
  en = energy_partition(x, v, t, box, expansion, index)
  for T in np.linspace(3.0, 0.0, 10001)[:-1]:
    for i in [3, 2, 1]:
      index_new = perturbate_system(index, i)
      en_new = energy_partition(x, v, t, box, expansion, index_new)
      de = en_new - en
      if de < 0:
        index = index_new
        en = en_new
      elif random.random() < np.exp(-de/T):
        index = index_new
        en = en_new
      print T, en
  value = np.zeros((npart + 1, 3))
  value[:, 0] = range(npart + 1)
  ec = index.copy()
  for clus in np.unique(ec):
    mass = np.count_nonzero(ec == clus)
    protons = np.count_nonzero((ec == clus) & (t.flatten() == 2))
    frac = float(protons)/mass
    value[mass, 1] += 1
    value[mass, 2] *= (value[mass, 1] - 1.0)/value[mass, 1]
    value[mass, 2] += frac/value[mass, 1]

  return value, (ec, [])

def brute_force(x, v, t, box, expansion):
  """
  Iterate through all possible partitions of the set and find the
  minimum

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

  expansion : float, optional
      The expansion velocity of the walls of the box, in box units.

  Returns
  -------

  Returns
  -------

  value, (ec, inf) : numpy array, numpy array, list
      value is the [mass, occupancy, fraction] histogram
      mst is the array of indices to which each particle belongs.
      inf is the list of infinite clusters, which is always empty.
  """
  npart = np.shape(x)[0]
  min_idx = np.arange(npart)
  min_en = energy_partition(x, v, t, box, expansion, min_idx)
  for part in partition(range(npart)):
    idx = np.arange(npart)
    for i, dx in enumerate(part):
      for j in dx:
        idx[j] = i
    en = energy_partition(x, v, t, box, expansion, idx)
    if en < min_en:
      min_en = en
      min_idx = idx

  value = np.zeros((npart + 1, 3))
  value[:, 0] = range(npart + 1)
  ec = min_idx.copy()
  for clus in np.unique(ec):
    mass = np.count_nonzero(ec == clus)
    protons = np.count_nonzero((ec == clus) & (t.flatten() == 2))
    frac = float(protons)/mass
    value[mass, 1] += 1
    value[mass, 2] *= (value[mass, 1] - 1.0)/value[mass, 1]
    value[mass, 2] += frac/value[mass, 1]

  return value, (ec, [])

def energy_partition(x, v, t, box, expansion, idx):
  """
  Get the total energy of a specific cluster partition

  Parameters
  ----------

  x : numpy float64 array
      Positions of the particles in the system

  v : numpy float64 array
      Velocities of the particles

  t : numpy int32 array
      Types of the particles

  idx : iterable
      Cluster index of each particle

  box : numpy float64 array
      Box

  expansion : float, optional
      The expansion velocity of the walls of the box, in box units.

  Returns
  -------

  energy : float
      Energy of the partition
  """
  size_x = box[0][1] - box[0][0]
  size_y = box[1][1] - box[1][0]
  size_z = box[2][1] - box[2][0]
  if size_x != size_y or size_y != size_z:
    raise ValueError("The box should be cubic for this to work")
  else:
    size = size_x

  natoms = np.shape(x)[0]
  x_p = x.ctypes.data_as(ct.c_void_p)
  v_p = v.ctypes.data_as(ct.c_void_p)
  t_p = t.ctypes.data_as(ct.c_void_p)
  idx_p = idx.ctypes.data_as(ct.c_void_p)
  enpart_c.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p,
                       ct.c_int, ct.c_double, ct.c_double,
                       ct.c_void_p]
  enpart_c.restype = ct.c_double;
  return enpart_c(x_p, v_p, t_p, natoms, size, expansion, idx_p)

def perturbate_system(idx, nperm=1):
  """
  Perform a perturbation of the cluster partition of the system

  Parameters
  ----------

  idx : iterable
      Cluster index of each particle

  nperm : integer, optional
      Number of cluster permutations to perform

  Returns
  -------

  idx_new : iterable
      Perturbation of the system
  """
  idx_new = copy.deepcopy(idx)
  pos = range(len(idx))
  orig = np.random.choice(pos, nperm)
  avail = np.unique(idx)
  empty = np.setdiff1d(pos, avail, assume_unique=True)
  if len(empty) != 0:
    avail = np.hstack((avail, empty))
  new = np.random.choice(avail, nperm)
  new = np.random.choice(pos, nperm)
  for o, n in zip(orig, new):
    idx_new[o] = n
  return idx_new
