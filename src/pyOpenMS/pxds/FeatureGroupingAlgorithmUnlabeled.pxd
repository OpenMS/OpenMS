from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from FeatureGroupingAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmUnlabeled(FeatureGroupingAlgorithm) :
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        FeatureGroupingAlgorithmUnlabeled() nogil except +
        # private
        FeatureGroupingAlgorithmUnlabeled(FeatureGroupingAlgorithmUnlabeled &) nogil except + # wrap-ignore
        void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out) nogil except +
        # POINTER # FeatureGroupingAlgorithm * create() nogil except +
        String getProductName() nogil except +
        void addToGroup(int map_id, FeatureMap feature_map) nogil except +
        void setReference(int map_id, FeatureMap map) nogil except +
        ConsensusMap getResultMap() nogil except +

