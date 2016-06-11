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

_DIRNAME = os.path.dirname(__file__)
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

def ecra(x, v, t):
  """
  Calculate cluster recognition algorithm with Simulated Annealing as
  in [Dorso]_

  .. [Dorso] Dorso and Randrup, Phys. Lett. B, 301, 4, 328-333

  Parameters
  ----------

  x : numpy float64 array
      Positions of the particles in the system

  v : numpy float64 array
      Velocities of the particles

  t : numpy int32 array
      Types of the particles

  Returns
  -------

  idx : iterable
      Cluster indices of the particles

  """
  npart = np.shape(x)[0]
  idx = range(npart)
  en = energy_partition(x, v, t, idx)
  for T in np.linspace(3.0, 0.0, 1001)[:-1]:
    idx_new = perturbate_system(idx)
    en_new = energy_partition(x, v, t, idx_new)
    de = en_new - en
    print idx, idx_new, de, en
    if de < 0:
      idx = idx_new
      en = en_new
    elif random.random() < np.exp(-de/T):
      idx = idx_new
      en = en_new
  return idx, en

def brute_force(x, v, t):
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

  Returns
  -------

  idx : iterable
      Cluster indices of the particles
  """
  npart = np.shape(x)[0]
  min_idx = range(npart)
  min_en = energy_partition(x, v, t, min_idx)
  for part in partition(range(npart)):
    idx = range(npart)
    for i, dx in enumerate(part):
      for j in dx:
        idx[j] = i
    en = energy_partition(x, v, t, idx)
    if en < min_en:
      min_en = en
      min_idx = idx
  return min_idx, min_en

def energy_partition(x, v, t, idx):
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

  Returns
  -------

  energy : float
      Energy of the partition
  """
  e = 0
  for clus in np.unique(idx):
    x_clus = x[idx == clus]
    v_clus = v[idx == clus]
    t_clus = t[idx == clus]
    nclus = np.sum(idx == clus)
    vcm = np.zeros(3)
    for i in range(nclus):
      vcm += v_clus[i]
      for j in range(i + 1, nclus):
        dx = np.linalg.norm(x_clus[i] - x_clus[j])
        e += potential(dx, t_clus[i], t_clus[j])
    vcm /= nclus
    for i in range(nclus):
      k = np.dot((v_clus[i] - vcm), (v_clus[i] - vcm)) * 938.0 / 2
      e += k
  return e

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
