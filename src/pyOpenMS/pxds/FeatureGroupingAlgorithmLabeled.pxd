from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureGroupingAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmLabeled(FeatureGroupingAlgorithm) :
        # wrap-doc:
        #  A map feature grouping algorithm for labeling techniques with two
        #  labels
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        FeatureGroupingAlgorithmLabeled() except + nogil 
        # private
        FeatureGroupingAlgorithmLabeled(FeatureGroupingAlgorithmLabeled &) except + nogil  # wrap-ignore
        void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out) except + nogil 
        # POINTER # FeatureGroupingAlgorithm * create() except + nogil 


