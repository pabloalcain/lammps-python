import numpy as np
import itertools as it
from postprocess.extract import Extraction


def replicate(x, t, box, nrep=1):
  """
  Replicate positions of particles in a box according to periodic
  boundary conditions.

  Parameters
  ----------

  x : numpy array
      Positions of the particles

  t : numpy array
      Types of the particles

  box : iterable
      Simulation box

  nrep : int, optional
      Number of repetitions to consider. Default value is 1.

  Returns
  -------

  x, t, box : numpy array, numpy array, box
      Positions, type and box of the repeated values.
  """
  boxnew = box.copy()
  rep = np.linspace(0, 2*nrep, 2*nrep + 1) - nrep
  for (orig, new) in zip(box, boxnew):
    new[0] -= nrep * (orig[1] - orig[0])
    new[1] += nrep * (orig[1] - orig[0])
  szorig = np.zeros(3)
  for i, _ in enumerate(szorig):
    szorig[i] = box[i][1] - box[i][0]
  norig = np.shape(x)[0]
  nnew = norig * (2*nrep + 1)**3
  xnew = np.zeros((nnew, 3))
  tnew = np.zeros((nnew, 1), dtype=np.int32)

  idx = 0
  for n, r in enumerate(x):
    for nx, ny, nz in it.product(rep, rep, rep):
      desp = np.array((nx, ny, nz)) * szorig
      xnew[idx] = np.array(r) + desp
      tnew[idx] = t[n]
      idx += 1
  return xnew, tnew, boxnew
