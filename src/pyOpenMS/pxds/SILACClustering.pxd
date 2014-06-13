from Types cimport *
from HierarchicalClustering cimport *
from SILACPattern cimport *

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/SILACClustering.h>" namespace "OpenMS":
    
    cdef cppclass SILACClustering:
        # inherits:
        #  HierarchicalClustering< SILACPattern * >
        SILACClustering(SILACClustering) nogil except + #wrap-ignore
        double rt_min # wrap-ignore
        double rt_max_spacing # wrap-ignore
        # SILACClustering(PointCoordinate & cluster_dimension, double rt_min, double rt_max_spacing) nogil except +
        void cluster() nogil except +

