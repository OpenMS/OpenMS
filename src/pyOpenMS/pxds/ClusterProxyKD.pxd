from Types cimport *
from libcpp cimport bool
from FeatureGroupingAlgorithm cimport *
from KDTreeFeatureMaps cimport *
from FeatureDistance cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>" namespace "OpenMS":
    
    cdef cppclass ClusterProxyKD "OpenMS::ClusterProxyKD":
        ClusterProxyKD() except + nogil 
        ClusterProxyKD(ClusterProxyKD &) except + nogil 
        ClusterProxyKD(Size size, double avg_distance, Size center_index) except + nogil 
        bool operator<(ClusterProxyKD & rhs) except + nogil 
        bool operator!=(ClusterProxyKD & rhs) except + nogil 
        bool operator==(ClusterProxyKD & rhs) except + nogil 
        Size getSize() except + nogil 
        bool isValid() except + nogil 
        double getAvgDistance() except + nogil 
        Size getCenterIndex() except + nogil 

