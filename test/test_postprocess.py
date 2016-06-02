from postprocess.extract import Extraction
from postprocess import tools
import numpy as np
import nose.tools as nst
import itertools as it

class TestExtract(object):
  """
  Test methods in the extract
  """

  def setUp(self):
    """Setup method with a mock extraction """
    self.ext = Extraction('data/')

  def test_init(self):
    """Clean init"""
    e = Extraction('data/')
    nst.assert_equal(e.path, 'data/')

  def test_init_no_file(self):
    """Init with no folder"""
    nst.assert_raises(IOError, Extraction, "nodata/")

  def test_return_entry_onetemp(self):
    """One temperature returns the only entry that matches"""
    entries = self.ext.entries({'T': 0.5})
    nst.assert_equal([_['id'] for _ in entries],
                     ['l20/N50/exp0.0/medium/x0.5/d0.05/T0.5/'])

  def test_return_entry_nomatch(self):
    """No match returns returns empty list"""
    entries = self.ext.entries({'x': 0.3})
    nst.assert_equal(entries, [])

  def test_return_entry_multiple_match(self):
    """Multiple match return whole list"""
    entries = self.ext.entries({'x': 0.5})
    nst.assert_equal(len(entries), 31)

  def test_particle_extraction(self):
    """Extract columns from an entry"""
    inf = self.ext.particle((2, 3, 4), {'T': 0.5})
    known = [-0.396773, -1.25902, -1.54895]
    for i, j in zip(inf[0][0, :], known):
      nst.assert_almost_equal(i, j)

  def test_mult_particle_extraction(self):
    """Extract columns from an entry with multiple hits"""
    inf = self.ext.particle((2, 3, 4), {'x': 0.5})
    nst.assert_equal(np.shape(inf), (31, 50, 3))

  def test_position_extraction(self):
    """Extract positions from an entry"""
    inf = self.ext.x({'T': 0.5}, idx=1)
    known = [-0.43157, -0.990112, -1.36708]
    for i, j in zip(inf[0][0, :], known):
      nst.assert_almost_equal(i, j)

  def test_velocity_extraction(self):
    """Extract velocities from an entry"""
    inf = self.ext.v({'T': 0.5}, idx=3)
    known = [-0.0167678, -0.0189928, 0.0155681]
    for i, j in zip(inf[0][0, :], known):
      nst.assert_almost_equal(i, j)

  def test_type_extraction(self):
    """Extract types from an entry"""
    inf = self.ext.t({'T': 0.5}, idx=3)
    known = [1]*25 + [2]*25
    for i, j in zip(inf[0][0, :], known):
      nst.assert_equal(i, j)

  def test_box_extraction(self):
    """Extract box from an entry"""
    box = self.ext.box({'T': 0.5}, idx=3)[0]
    known = [[-5, 5], [-5, 5], [-5, 5]]
    for i, j in zip(box, known):
      nst.assert_equal(i[0], j[0])
      nst.assert_equal(i[1], j[1])

class TestTools(object):
  """
  Class to test the tools.
  """

  def test_many_part_replication(self):
    """
    Replication returns expected size
    """
    nrep = 1
    x = np.array([[1, 2, 3], [4, 5, 1], [2, 4, 3], [0.5, 2, -0.5],
                  [1.8, 6, 1.0], [-2, 5, -1]])
    t = np.array([[1], [2], [2], [2], [2], [2]], dtype=np.int32)
    box = np.array([[0, 10], [0, 6], [-1, 4]])
    boxrep = np.array([[-10, 20], [-6, 12], [-6, 9]])
    xn, tn, boxn = tools.replicate(x, t, box, nrep=nrep)
    nst.assert_equal(len(xn), len(x)*(2*nrep+1)**3)
    nst.assert_equal(len(tn), len(t)*(2*nrep+1)**3)
    for i, j in zip(boxn, boxrep):
      nst.assert_equal(i[0], j[0])
      nst.assert_equal(i[1], j[1])

  def test_one_part_replication(self):
    """
    Compare replication with hand-made replicas
    """
    nrep = 1
    x = np.array([[0.0, 0.0, 0.0]])
    t = np.array([[1]], dtype=np.int32)
    trep = np.array([[1]]*27, dtype=np.int32)
    box = np.array([[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]])
    xn, tn, boxn = tools.replicate(x, t, box, nrep=nrep)
    nst.assert_equal(len(tn), len(trep))
    for i, j in zip(trep, tn):
      nst.assert_equal(i, j)
    rep = np.array((-1, 0, 1))
    idx = 0
    for b in it.product(rep, rep, rep):
      for i in range(3):
        nst.assert_equal(xn[idx, i], b[i])
      idx += 1
