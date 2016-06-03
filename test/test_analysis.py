import analysis as A
import numpy as np
from postprocess import tools as T
from postprocess.extract import Extraction
import nose.tools as nst
import warnings

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
    #nst.assert_almost_equal(sum(abs(sk_rep.flatten() - sk.flatten())), 0)


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
    box = e.box({'T': 2.0})[0]
    size = box[0][1] - box[0][0]
    pairs = (((0,), (0,),),)
    gr = A.rdf(x, t, size, pairs, 100, False)
    d = np.shape(x)[0]/(float(size) ** 3)
    sk_transf = A.RDF.RDF.ssf(gr, d, False)
    k = sk_transf[:, 0].copy()
    sk = A.structureFactor(x, t, size, pairs, k, 1, lebedev=194)
    #nst.assert_almost_equal(sum(abs(sk_transf.flatten() - sk.flatten())), 0)

class TestMSTE(object):
  """
  Test Minimum Spanning Tree in analysis module
  """

  def test_mste_exists(self):
    """
    Can call mste
    """
    A.mste

  def test_graph_partition(self):
    """
    Returns partition of a graph.
    """
    graph = {1: [2, 3], 2: [1,], 3: [1,],
             4: [5, 6], 5: [6, 4], 6: [4, 5]}
    hand_part = [set([1, 2, 3]), set([4, 5, 6])]
    part = A.MSTE.MSTE._partition(graph)
    for i, j in zip(part, hand_part):
      nst.assert_equal(i, j)
    graph = {1: [2, 3], 2: [1,], 3: [1,],
             4: [5], 5: [4], 6: []}
    hand_part = [set([1, 2, 3]), set([4, 5]), set([6,])]
    part = A.MSTE.MSTE._partition(graph)
    for i, j in zip(part, hand_part):
      nst.assert_equal(i, j)

  def test_idx(self):
    """
    Binary to wall conversion and back
    """
    nst.assert_equal(A.MSTE.MSTE._idx2wall(1), "W")
    nst.assert_equal(A.MSTE.MSTE._idx2wall(-1), "E")
    nst.assert_equal(A.MSTE.MSTE._idx2wall(2), "N")
    nst.assert_equal(A.MSTE.MSTE._idx2wall(-2), "S")
    nst.assert_equal(A.MSTE.MSTE._idx2wall(4), "U")
    nst.assert_equal(A.MSTE.MSTE._idx2wall(-4), "D")
    nst.assert_equal(1, A.MSTE.MSTE._wall2idx("W"))
    nst.assert_equal(-1, A.MSTE.MSTE._wall2idx("E"))
    nst.assert_equal(2, A.MSTE.MSTE._wall2idx("N"))
    nst.assert_equal(-2, A.MSTE.MSTE._wall2idx("S"))
    nst.assert_equal(4, A.MSTE.MSTE._wall2idx("U"))
    nst.assert_equal(-4, A.MSTE.MSTE._wall2idx("D"))
    nst.assert_raises(ValueError, A.MSTE.MSTE._wall2idx, "GF")
    nst.assert_raises(ValueError, A.MSTE.MSTE._idx2wall, -8)


  def test_graph_creation(self):
    """
    Create a graph structure from connections
    """
    conn = [[1, 2, 1], [2, 4, 2], [1, 2, -2], [3, 3, -1]]
    hand_graph = {1: [2,], 2: [4,], 3: [3,]}
    hand_cnct = {(1, 2): ['W', 'S'], (2, 4): ['N'], (3, 3): ['E']}
    graph, cnct = A.MSTE.MSTE._create_graph(conn)
    nst.assert_equal(hand_graph, graph)
    nst.assert_equal(hand_cnct, cnct)

  def test_cycles_finding(self):
    """
    Return cycles and edges in undirected graph
    """
    graph = {1: [2, 4], 2: [1, 4,], 3: [3,], 4: [1, 2, 5], 5: [4,]}
    cycles = A.MSTE.MSTE._find_cycles(graph)
    hand_cycles = [[1, 2], [1, 4], [2, 4], [4, 5], [3,], [2, 4, 1]]
    for i in cycles:
      nst.assert_in(i, hand_cycles)

  def test_path(self):
    """
    Paths in the graph and whether they are infinite
    """
    #Connections here are duplicated because it's easier to implement
    #the c-code that generates them
    conn = [[1, 2, 1], [2, 1, -1], [2, 4, 2], [4, 2, -2], [4, 1, -1], [1, 4, 1],
            [4, 5, -4], [5, 4, 4], [3, 3, 1]]
    graph, cnct = A.MSTE.MSTE._create_graph(conn)
    paths, inf = A.MSTE.MSTE._find_paths(graph, cnct)
    hand_paths = [([2, 4, 1, 2], 'NEW'), ([3, 3], 'W'), ([1, 2, 1], 'WE'),
                  ([1, 4, 1], 'WE'), ([2, 4, 2], 'NS'), ([4, 5, 4], 'DU')]
    hand_inf = [[2, 4, 1, 2], [3, 3]]
    nst.assert_equal(hand_paths, paths)
    nst.assert_equal(hand_inf, inf)


  def test_cluster_1(self):
    """
    MST Cluster recognition: 1 cluster [no energy considerations]
    """
    x = np.array([[1, 2, 3], [4, 5, 1], [2, 4, 3], [0.5, 2, -0.5],
                  [1.8, 5, 1.0], [2, 5, 1]])
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 100
    index = A.MSTE.MSTE.cluster(x, v, t, size, False)
    hand_index = [0, 0, 0, 0, 0, 0]
    for i, j in zip(index, hand_index):
      nst.assert_equal(i, j)

  def test_cluster_2(self):
    """
    MST Cluster recognition: 2 clusters [no energy considerations]
    """

    x = np.array([[1.8, 55, 3], [4, 5, 1], [2, 4, 3], [0.5, 2, -0.5],
                  [1.8, 56, 1.0], [2, 5, 1]])
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 100
    index = A.MSTE.MSTE.cluster(x, v, t, size, False)
    hand_index = [0, 1, 1, 1, 0, 1]
    for i, j in zip(index, hand_index):
      nst.assert_equal(i, j)

  def test_cluster_3(self):
    """
    MST Cluster recognition: 6 clusters [no energy considerations]
    """

    x = np.array([[1.8, 55, 3], [54, 55, 1], [2, 54, 53], [0.5, 25, 0.5],
                  [1.8, 26, 51.0], [22, 25, 1]])
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 100
    index = A.MSTE.MSTE.cluster(x, v, t, size, False)
    hand_index = [0, 1, 2, 3, 4, 5]
    for i, j in zip(index, hand_index):
      nst.assert_equal(i, j)

  def test_connections_1(self):
    """
    Connections: No connections in the cluster
    """
    x = np.array([[1.8, 55, 3], [54, 55, 1], [2, 54, 53], [0.5, 25, 0.5],
                  [1.8, 26, 51.0], [22, 25, 1]])
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 100
    index = A.MSTE.MSTE.cluster(x, v, t, size, False)
    conn = A.MSTE.MSTE.connections(x, v, t, index, size, False, 0.0)
    nst.assert_equal(conn.shape, (0, 3))


  def test_connections_2(self):
    """
    Connections: Cluster 0 connects with 4
    """
    x = np.array([[1.8, 55, 3], [54, 55, 1], [2, 54, 53], [0.5, 25, 0.5],
                  [99.8, 55, 3], [22, 25, 1]])
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 100.0
    index = A.MSTE.MSTE.cluster(x, v, t, size, False)
    conn = A.MSTE.MSTE.connections(x, v, t, index, size, False, 0.0)
    nst.assert_equal(conn.shape, (2, 3))

  def test_connections_3(self):
    """
    Connections: Cluster 0 connects twice with 4
    """
    x = np.array([[1.8, 0.5, 3], [54, 55, 1], [2, 54, 53], [0.5, 25, 0.5],
                  [99.8, 99.5, 3], [22, 25, 1]])
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 100.0
    index = A.MSTE.MSTE.cluster(x, v, t, size, False)
    conn = A.MSTE.MSTE.connections(x, v, t, index, size, False, 0.0)
    nst.assert_equal(conn.shape, (4, 3))

  def test_connections_4(self):
    """
    Connections: Cluster 0 connects with itself
    """
    x = np.array([[1, 2, 3], [2, 2, 3], [3, 2, 3], [5, 2, 3],
                  [7, 2, 3], [9, 2, 3]], dtype=np.float64)
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 10
    index = A.MSTE.MSTE.cluster(x, v, t, size, False)
    conn = A.MSTE.MSTE.connections(x, v, t, index, size, False, 0.0)
    nst.assert_equal(conn.shape, (2, 3))

  def test_mst_1(self):
    """
    Full MST Cluster warns when x, v and t are not the proper type
    """
    x = np.array([[1, 2, 3], [2, 2, 3], [3, 2, 3], [5, 2, 3],
                  [7, 2, 3], [9, 2, 3]], dtype=np.int32)
    v = np.zeros((6, 3), dtype=np.int32)
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.float64)
    size = 10
    with warnings.catch_warnings(record=True) as w:
      value, _ = A.mste(x, v, t, size, False)
    nst.assert_equal(len(w), 3)

  def test_mst_2(self):
    """
    Full MST Cluster with tallies: 1 cluster [no energy considerations]
    """
    x = np.array([[1, 2, 3], [4, 5, 1], [2, 4, 3], [0.5, 2, 0.5],
                  [1.8, 5, 1.0], [2, 5, 1]])
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 100
    value, _ = A.mste(x, v, t, size, False)
    value_hand = np.zeros((7, 3))
    value_hand[:, 0] = range(7)
    value_hand[6, 1:] = (1, 0.5)
    for i, j in zip(value, value_hand):
      for e1, e2 in zip(i, j):
        nst.assert_equal(e1, e2)
    print value

  def test_mst_3(self):
    """
    Full MST Cluster with tallies: 1 cluster links with itself
    """
    x = np.array([[1, 2, 3], [2, 2, 3], [3, 2, 3], [5, 2, 3],
                  [7, 2, 3], [9, 2, 3]], dtype=np.float64)
    v = np.zeros((6, 3))
    t = np.array([[1], [2], [1], [2], [1], [2]], dtype=np.int32)
    size = 10
    value, _ = A.mste(x, v, t, size, False)
    value_hand = np.zeros((7, 3))
    value_hand[:, 0] = range(7)
    value_hand[0, 1:] = (1, 0.5)
    for i, j in zip(value, value_hand):
      for e1, e2 in zip(i, j):
        nst.assert_equal(e1, e2)

  def test_mst_4(self):
    """
    Full MST Cluster error when cast of t modifies its value
    """
    x = np.array([[1, 2, 3], [2, 2, 3], [3, 2, 3], [5, 2, 3],
                  [7, 2, 3], [9, 2, 3]], dtype=np.int32)
    v = np.zeros((6, 3), dtype=np.int32)
    t = np.array([[1.5], [2], [1], [2], [1], [2]], dtype=np.float64)
    size = 10
    nst.assert_raises(ValueError, A.mste, x, v, t, size, False)
