from DefaultParamHandler cimport *
from ConsensusMap cimport *
from FeatureMap cimport *
from Feature cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FeatureGroupingAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler
        # wrap-ignore
        void transferSubelements(libcpp_vector[ConsensusMap] maps,
                                 ConsensusMap & out
                                ) nogil except +

        # since this is a base class, cannot have overloaded methods
        # void group(libcpp_vector[ FeatureMap[Feature] ] & maps, ConsensusMap & out)
        # void group(libcpp_vector[ ConsensusMap ] & maps, ConsensusMap & out)
        void registerChildren()
