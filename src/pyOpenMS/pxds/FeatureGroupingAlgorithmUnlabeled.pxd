from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from FeatureGroupingAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmUnlabeled(FeatureGroupingAlgorithm) :
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        FeatureGroupingAlgorithmUnlabeled() except + nogil 
        # private
        FeatureGroupingAlgorithmUnlabeled(FeatureGroupingAlgorithmUnlabeled &) except + nogil  # wrap-ignore
        void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out) except + nogil 
        # POINTER # FeatureGroupingAlgorithm * create() except + nogil 
        void addToGroup(int map_id, FeatureMap feature_map) except + nogil 
        void setReference(int map_id, FeatureMap map) except + nogil 
        ConsensusMap getResultMap() except + nogil 

