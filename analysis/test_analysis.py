"""
Tests for analysis
"""

import unittest
import structure
import numpy as np

def fun(x):
  return x + 1

class UnitTesting(unittest.TestCase):
  def setUp(self):
    x = np.array([[1, 1, 1], [2, 1, 1], [3, 1, 1],
                  [4, 1, 1], [5, 1, 1], [6, 1, 1]], dtype=np.double)
    t = np.array([1, 1, 2, 2, 1, 2], dtype=np.int32)
    size = 10
    self.g_test = structure.rdf(x, t, size)

  def test_rdf_default(self):
    x = np.array([[1, 1, 1], [2, 1, 1], [3, 1, 1],
                  [4, 1, 1], [5, 1, 1], [6, 1, 1]], dtype=np.double)
    t = np.array([1, 1, 2, 2, 1, 2], dtype=np.int32)
    size = 10
    g = structure.rdf(x, t, size)
    self.assertEqual(np.shape(g), (200, 5))

  def test_rdf_casting(self):
    x = np.array([[1, 1, 1], [2, 1, 1], [3, 1, 1],
                  [4, 1, 1], [5, 1, 1], [6, 1, 1]], dtype=np.int)
    t = np.array([1, 1, 2, 2, 1, 2], dtype=np.double)
    size = 10
    g_cast = structure.rdf(x, t, size)
    self.assertTrue((g_cast == self.g_test).all())

  def test_ssf(self):
    density = 0.006
    s = structure.ssf(self.g_test[:, [0, 1]], density, npoints=200)
    self.assertEqual(np.shape(s[0]), (200, 2))

  def test_ssf_def(self):
    x = np.array([[1, 1, 1], [2, 1, 1], [3, 1, 1],
                  [4, 1, 1], [5, 1, 1], [6, 1, 1]], dtype=np.int)
    t = np.array([1, 1, 2, 2, 1, 2], dtype=np.double)
    size = 10
    k = np.linspace(0, 2*np.pi*200/size, 200)
    density = 0.006
    s = structure.ssf(self.g_test[:, [0, 1]], density, npoints=200)
    k = s[0][:, 0].copy()
    ssf = structure.ssf_def(x, t, k, size)
    import pylab
    pylab.plot(ssf[:, 0], ssf[:, 1])
    pylab.plot(s[0][:, 0], s[0][:, 1])
    pylab.title("These two should be equal")
    pylab.show()


if __name__ == '__main__':
  unittest.main()
