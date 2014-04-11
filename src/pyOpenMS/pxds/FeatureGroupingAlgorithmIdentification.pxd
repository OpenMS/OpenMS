from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureGroupingAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmIdentification.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmIdentification(FeatureGroupingAlgorithm) :
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        FeatureGroupingAlgorithmIdentification() nogil except +
        FeatureGroupingAlgorithmIdentification(FeatureGroupingAlgorithmIdentification) nogil except + #wrap-ignore
        void group(libcpp_vector[ FeatureMap[Feature] ] & maps, ConsensusMap & out) nogil except +
        # POINTER # FeatureGroupingAlgorithm * create() nogil except +
        String getProductName() nogil except +

