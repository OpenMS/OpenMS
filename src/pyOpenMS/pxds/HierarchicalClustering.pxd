from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/HierarchicalClustering.h>" namespace "OpenMS":
    
    cdef cppclass HierarchicalClustering[PointRef]:
        HierarchicalClustering() nogil except + #wrap-ignore
        HierarchicalClustering(HierarchicalClustering) nogil except + #wrap-ignore
        # Grid grid
        HierarchicalClustering(DPosition2 & cluster_dimension) nogil except +
        # NAMESPACE # Grid::cell_iterator insertPoint(PointCoordinate & d, PointRef & ref) nogil except +
        void cluster() nogil except +

