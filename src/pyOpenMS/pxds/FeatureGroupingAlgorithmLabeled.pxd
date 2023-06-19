from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureGroupingAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmLabeled(FeatureGroupingAlgorithm) :
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        FeatureGroupingAlgorithmLabeled() nogil except +
        # private
        FeatureGroupingAlgorithmLabeled(FeatureGroupingAlgorithmLabeled &) nogil except + # wrap-ignore
        void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out) nogil except +
        # POINTER # FeatureGroupingAlgorithm * create() nogil except +
        String getProductName() nogil except +

