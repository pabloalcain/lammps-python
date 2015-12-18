import extract as E
import tools as T
import numpy as np
import analysis as A
import pylab as pl

def test_replicas_scheme(full=False):
  """In this test, we check that the Replica scheme we use is equal to
  the result with the replicas made ad-hoc in a new file"""

  k1 = np.linspace(0.0, 0.5, 100)
  k2 = np.linspace(2.0, 4.0, 100)
  k = np.hstack((k1, k2))
  T.replicate('test.lammpstrj', 'replica.lammpstrj')

  x = E.positions('test.lammpstrj', 0)
  t = E.types('test.lammpstrj', 0)
  b = E.box('test.lammpstrj')
  sz = b[0][1] - b[0][0]
  sk = A.ssf(x, t, k, sz, 1)

  x = E.positions('replica.lammpstrj', 0)
  t = E.types('replica.lammpstrj', 0)
  b = E.box('replica.lammpstrj')
  sz = b[0][1] - b[0][0]
  sk_rep = A.ssf(x, t, k, sz, 0)

  print 'Difference should be close to the float eps: sum(abs(sk_rep - sk)) < 1e-8*len(k)'
  print 'sum(abs(sk_rep - sk)) = {0}' .format(sum(abs(sk_rep - sk)))

  if full:
    pl.plot(sk[:, 0], sk[:, 1], label='Replica collapse')
    pl.plot(sk_rep[:, 0], sk_rep[:, 1], label='Replicas ad-hoc')
    pl.legend()
    pl.title('For the test to pass, both should be the same graph')
    pl.show()
if __name__ == '__main__':
  test_replicas_scheme(True)
