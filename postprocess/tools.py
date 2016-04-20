import numpy as np
import itertools as it
import extract

def replicate(fname, fout, nrep=1):
  pos = extract.positions(fname, 0)
  typ = extract.types(fname, 0)
  boxorig = extract.box(fname)
  boxnew = boxorig.copy()

  rep = np.linspace(0, 2*nrep, 2*nrep + 1) - nrep
  for (orig, new) in zip(boxorig, boxnew):
    new[0] -= nrep * (orig[1] - orig[0])
    new[1] += nrep * (orig[1] - orig[0])
  szorig = np.zeros(3)
  for i, sz in enumerate(szorig):
    szorig[i] = boxorig[i][1] - boxorig[i][0]
  norig = np.shape(pos)[0]
  nnew = norig * (2*nrep + 1)**3
  f = open(fout, 'w')
  print>>f, "ITEM: TIMESTEP"
  print>>f, 0
  print>>f, "ITEM: NUMBER OF ATOMS"
  print>>f, nnew
  print>>f, "ITEM: BOX BOUNDS pp pp pp"
  for bx in boxnew:
    print>>f, bx[0], bx[1]
  print>>f, "ITEM: ATOMS id type x y z vx vy vz"

  idx = 0;
  for n, r in enumerate(pos):
    for nx, ny, nz in it.product(rep, rep, rep):
      idx += 1
      desp = np.array((nx, ny, nz)) * szorig
      print>>f, idx, typ[n],
      print>>f,  ' '.join(map(str, np.array(r) + desp)), '0.0 0.0 0.0'

if __name__ == "__main__":
  replicate("test.lammpstrj", "replica.lammpstrj")
