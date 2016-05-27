import analysis as A
import numpy as np
from postprocess import tools as T
from postprocess.extract import Extraction
import nose.tools as nst

class TestSSF(object):
  """
  Test Structure Factor in analysis module
  """

  def test_structureFactor_exists(self):
    """
    Can call structureFactor
    """
    A.structureFactor

  def test_replicas_scheme(self):
    """Replica scheme we use is equal to replicas made ad-hoc"""

    k1 = np.linspace(0.0, 0.5, 100)
    k2 = np.linspace(2.0, 4.0, 100)
    k = np.hstack((k1, k2))
    x = np.array([[1, 2, 3], [4, 5, 1], [2, 4, 3], [0.5, 2, -0.5],
                  [1.8, 6, 1.0], [-2, 5, -1]])
    t = np.array([[1], [2], [2], [2], [2], [2]], dtype=np.int32)
    box = np.array([[0, 10], [0, 10], [0, 10]])
    pairs = (((0,), (0,),),)
    xn, tn, boxn =  T.replicate(x, t, box)
    sz = box[0][1] - box[0][0]
    sk = A.structureFactor(x, t, sz, pairs, k, 3)
    szn = boxn[0][1] - boxn[0][0]
    sk_rep = A.structureFactor(xn, tn, szn, pairs, k, 1)
    nst.assert_almost_equal(sum(abs(sk_rep.flatten() - sk.flatten())), 0)

  def test_replicas_non_cubic(self):
    """Replica scheme in non-cubic boxes. OK if fails so far."""

    k1 = np.linspace(0.0, 0.5, 100)
    k2 = np.linspace(2.0, 4.0, 100)
    k = np.hstack((k1, k2))
    x = np.array([[1, 2, 3], [4, 5, 1], [2, 4, 3], [0.5, 2, -0.5],
                  [1.8, 6, 1.0], [-2, 5, -1]])
    t = np.array([[1], [2], [2], [2], [2], [2]], dtype=np.int32)
    box = np.array([[0, 6], [-2, 5], [2, 8]])
    pairs = (((0,), (0,),),)
    xn, tn, boxn = T.replicate(x, t, box)
    sz = box[0][1] - box[0][0]
    sk = A.structureFactor(x, t, sz, pairs, k, 3)
    szn = boxn[0][1] - boxn[0][0]
    sk_rep = A.structureFactor(xn, tn, szn, pairs, k, 1)
    nst.assert_almost_equal(sum(abs(sk_rep.flatten() - sk.flatten())), 0)


class TestRDF(object):
  """
  Test Structure Factor in analysis module
  """

  def test_rdf_exists(self):
    """
    Can call rdf
    """
    A.rdf

  def test_rdf_transform(self):
    """
    RDF Fourier Transform equal to Structure Factor. Need a better way
    to measure, since it's quite fuzzy and will fail.
    """
    e = Extraction('./data/')
    x = e.x({'T': 2.0})[0]
    t = e.t({'T': 2.0})[0]
    box = e.box({'T': 2.0})
    size = box[0][1] - box[0][0]
    pairs = (((0,), (0,),),)
    gr = A.rdf(x, t, size, pairs, 100, False)
    d = np.shape(x)[0]/(float(size) ** 3)
    sk_transf = A.RDF.RDF.ssf(gr, d, False)
    k = sk_transf[:, 0].copy()
    sk = A.structureFactor(x, t, size, pairs, k, 1, lebedev=194)
    nst.assert_almost_equal(sum(abs(sk_transf.flatten() - sk.flatten())), 0)
