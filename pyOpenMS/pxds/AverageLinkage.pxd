from libcpp.vector cimport vector as libcpp_vector
from ProgressLogger cimport *
# from DistanceMatrix cimport *
# from ClusterFunctor cimport *

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>" namespace "OpenMS":
    
    cdef cppclass AverageLinkage "OpenMS::AverageLinkage":
        AverageLinkage() nogil except +
        AverageLinkage(AverageLinkage) nogil except +
        # void operator()(DistanceMatrix[ float ] &original_distance, libcpp_vector[ BinaryTreeNode ] &cluster_tree, float threshold) nogil except +
        # ClusterFunctor * create() nogil except +
        String getProductName() nogil except +

