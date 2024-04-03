from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from FeatureGroupingAlgorithm cimport *
from KDTreeFeatureMaps cimport *
from FeatureDistance cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>" namespace "OpenMS":
    
    cdef cppclass FeatureGroupingAlgorithmKD(FeatureGroupingAlgorithm,ProgressLogger) :
        # wrap-inherits:
        #  FeatureGroupingAlgorithm
        #  ProgressLogger
        FeatureGroupingAlgorithmKD() except + nogil  # wrap-doc:A feature grouping algorithm for unlabeled data
        # private
        FeatureGroupingAlgorithmKD(FeatureGroupingAlgorithmKD &) except + nogil  # wrap-ignore
        void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out) except + nogil 
        void group(libcpp_vector[ ConsensusMap ] & maps, ConsensusMap & out) except + nogil 
        # POINTER # FeatureGroupingAlgorithm * create() except + nogil 
      