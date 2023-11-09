from libcpp.vector cimport vector as libcpp_vector
from ProgressLogger cimport *
# from DistanceMatrix cimport *
# from ClusterFunctor cimport *

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>" namespace "OpenMS":
    
    cdef cppclass AverageLinkage "OpenMS::AverageLinkage":
        AverageLinkage() except + nogil 
        AverageLinkage(AverageLinkage &) except + nogil 
        # void operator()(DistanceMatrix[ float ] &original_distance, libcpp_vector[ BinaryTreeNode ] &cluster_tree, float threshold) except + nogil 
        # ClusterFunctor * create() except + nogil 
        String getProductName() except + nogil 

