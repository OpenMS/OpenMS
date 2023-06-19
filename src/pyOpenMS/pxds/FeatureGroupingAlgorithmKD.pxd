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
        FeatureGroupingAlgorithmKD() nogil except + # wrap-doc:A feature grouping algorithm for unlabeled data
        # private
        FeatureGroupingAlgorithmKD(FeatureGroupingAlgorithmKD &) nogil except + # wrap-ignore
        void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out) nogil except +
        void group(libcpp_vector[ ConsensusMap ] & maps, ConsensusMap & out) nogil except +
        # POINTER # FeatureGroupingAlgorithm * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Returns the product name (for the Factory)

