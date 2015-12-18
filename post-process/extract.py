"""
This file has routines that extract information from dump lammps files
"""
import numpy as np

def positions(fname, step):
  """
  Extract positions at a definite timestep from a lammps file
  """
  f = open(fname)
  n = 0
  for line in f.readlines():
    n += 1
    if n < 4:
      continue
    if n == 4:
      npart = int(line)
      s1 = 9 * (step + 1) +  step * npart + 1
      x = np.zeros((npart, 3))
      continue
    if n < s1:
      continue
    idxp = n - s1
    for j in range(3): x[idxp, j] = line.split()[j + 2]
    if n == s1 + npart -1:
      break
  f.close()
  return x

def types(fname, step):
  """
  Extract types at a definite timestep from a lammps file
  """
  f = open(fname)
  n = 0
  for line in f.readlines():
    n += 1
    if n < 4:
      continue
    if n == 4:
      npart = int(line)
      s1 = 9 * (step + 1) +  step * npart + 1
      t = np.zeros(npart, dtype=np.int32)
      continue
    if n < s1:
      continue
    idxp = n - s1
    t[idxp] = line.split()[1]
    if n == s1 + npart -1:
      break
  f.close()
  return t

def box(fname):
  """
  Extract box at the first timestep from a lammps file
  """
  f = open(fname)
  size = np.zeros((3, 2))
  n = 0
  for line in f.readlines():
    n += 1
    if n > 5 and n < 9:
      size[n-6, :] = line.split()
      continue
    if n == 9:
      break
  f.close()
  return size

if __name__ == "__main__":
  x = positions("dump.lammpstrj", 1)
  t = types("dump.lammpstrj", 1)
  size = box("test.lammpstrj")
