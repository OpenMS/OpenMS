from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureGroupingAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmLabeled(FeatureGroupingAlgorithm) :
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        FeatureGroupingAlgorithmLabeled() except + nogil 
        # private
        FeatureGroupingAlgorithmLabeled(FeatureGroupingAlgorithmLabeled &) except + nogil  # wrap-ignore
        void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out) except + nogil 
        # POINTER # FeatureGroupingAlgorithm * create() except + nogil 


