import analysis as A
import analysis.ECRA.ECRA as ECRA
from postprocess.extract import Extraction
import numpy as np
import nose.tools as nst
import warnings

class TestECRA(object):
  """
  Test Early Cluster Recognition module
  """

  def setUp(self):
    self.npart = 8
    self.x = np.array([[0, 0, 0.0], [0, 0, 1.4], [0, 0, 2.8],
                       [0, 0, 4.3], [0, 0, 5.6], [0, 0, 7.0],
                       [0, 0, 8.4], [0, 0, 9.8]],  dtype=np.float64)
    self.v = np.array([[0, 0, 0.000], [0, 0, 0.001], [0, 0, 0.002],
                       [0, 0, 0.003], [0, 0, 0.004], [0, 0, 0.005],
                       [0, 0, 0.006], [0, 0, 0.007]], dtype=np.float64)
    self.t = np.array([[1], [2], [1], [2], [1], [2], [1], [2]],
                      dtype=np.int32)

    self.v01 = -11.575673959850974
    self.v02 = 1.9982570114653591
    self.v23 = -11.422591404313323
    self.large_box = [[-100, 100], [-100, 100], [-100, 100]]
    self.t01 = 0.0002345
    self.t02 = 0.0009380
    self.t23 = 0.0002345

  def test_ECRA_exists(self):
    """
    Can call ecra
    """
    A.ecra

  def test_perturbation(self):
    """
    Can't do anything due to randomness.
    """
    idx = range(self.npart)
    idx_new = ECRA.perturbate_system(idx, 3)
    nst.assert_equal(len(idx), len(idx_new))

  def test_potential(self):
    """
    Test potential between two particles
    """
    v01 = ECRA.potential(1.4, 1, 2)
    v02 = ECRA.potential(2.8, 1, 1)
    v23 = ECRA.potential(1.5, 1, 2)
    nst.assert_almost_equal(v02, self.v02)
    nst.assert_almost_equal(v01, self.v01)
    nst.assert_almost_equal(v23, self.v23)

  def test_energy_partition_1(self):
    """
    Test that all particles in monoclusters are zero energy
    """
    idx = range(self.npart)
    e = ECRA.energy_partition(self.x, self.v, self.t, self.large_box, 0.0, idx)
    nst.assert_equal(e, 0)

  def test_energy_partition_2(self):
    """
    Test that the cluster of particle 0 and 1 has the correct energy.
    """
    idx = range(self.npart)
    idx[1] = 0
    e = ECRA.energy_partition(self.x, self.v, self.t, self.large_box, 0.0, idx)
    nst.assert_almost_equal(e, self.v01 + self.t01)

  def test_energy_partition_3(self):
    """
    Test that the cluster of particle 0 and 2 has the correct energy.
    """
    idx = range(self.npart)
    idx[2] = 0
    e = ECRA.energy_partition(self.x, self.v, self.t, self.large_box, 0.0, idx)
    nst.assert_almost_equal(e, self.v02 + self.t02)

  def test_energy_partition_4(self):
    """
    Test that two clusters are additive.
    """
    idx = range(self.npart)
    idx[1] = 0
    idx[3] = 2
    e = ECRA.energy_partition(self.x, self.v, self.t, self.large_box, 0.0, idx)
    nst.assert_almost_equal(e, self.v01 + self.v23 + self.t01 + self.t23)

  def test_ecra(self):
    """
    ECRA gives good results for silly toy model.
    """
    #idx, en = ECRA.ecra(self.x, self.v, self.t, self.large_box, 0.0)
    idx, en = ECRA.brute_force(self.x, self.v, self.t, self.large_box, 0.0)
    idx_hand = [0]*8
    for i, j in zip(idx, idx_hand):
      nst.assert_equal(i, j)
    #nst.assert_equal(0, (en, idx))
