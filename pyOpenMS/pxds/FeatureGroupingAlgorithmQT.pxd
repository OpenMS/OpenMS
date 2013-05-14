from FeatureGroupingAlgorithm cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Feature cimport *
from FeatureMap cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>" namespace "OpenMS":

    cdef cppclass FeatureGroupingAlgorithmQT(FeatureGroupingAlgorithm):
        # wrap-inherits:
        #    FeatureGroupingAlgorithm

        FeatureGroupingAlgorithmQT()
        # not in FeatureGroupingAlgorithm, as cython has a problem with
        # overloaded methods in base classes
        void group(libcpp_vector[FeatureMap[Feature]] & maps,
                   ConsensusMap & out
                   )

        void group(libcpp_vector[ConsensusMap] & maps,
                   ConsensusMap & out
                   )
