from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureGroupingAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmUnlabeled(FeatureGroupingAlgorithm) :
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        FeatureGroupingAlgorithmUnlabeled() nogil except +
        FeatureGroupingAlgorithmUnlabeled(FeatureGroupingAlgorithmUnlabeled) nogil except + #wrap-ignore
        void group(libcpp_vector[ FeatureMap[Feature] ] & maps, ConsensusMap & out) nogil except +
        # POINTER # FeatureGroupingAlgorithm * create() nogil except +
        String getProductName() nogil except +

