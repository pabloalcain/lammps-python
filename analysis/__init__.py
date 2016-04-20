#from analysis.cluster import cluster, connections, find_paths
#from analysis.structure import rdf, ssf

import ctypes as C

libanalysis = C.CDLL('libanalysis.so')
cluster_c = libanalysis.cluster
connections_c = libanalysis.connections
rdf_c = libanalysis.rdf
ssf_c = libanalysis.ssf
