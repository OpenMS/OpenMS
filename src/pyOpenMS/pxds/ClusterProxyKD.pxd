from Types cimport *
from libcpp cimport bool
from FeatureGroupingAlgorithm cimport *
from KDTreeFeatureMaps cimport *
from FeatureDistance cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>" namespace "OpenMS":
    
    cdef cppclass ClusterProxyKD "OpenMS::ClusterProxyKD":
        ClusterProxyKD() nogil except +
        ClusterProxyKD(ClusterProxyKD &) nogil except +
        ClusterProxyKD(Size size, double avg_distance, Size center_index) nogil except +
        bool operator<(ClusterProxyKD & rhs) nogil except +
        bool operator!=(ClusterProxyKD & rhs) nogil except +
        bool operator==(ClusterProxyKD & rhs) nogil except +
        Size getSize() nogil except +
        bool isValid() nogil except +
        double getAvgDistance() nogil except +
        Size getCenterIndex() nogil except +

