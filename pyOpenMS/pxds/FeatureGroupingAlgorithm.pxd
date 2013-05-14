from DefaultParamHandler cimport *
from ConsensusMap cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FeatureGroupingAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler
        # wrap-ignore
        void transferSubelements(libcpp_vector[ConsensusMap] maps,
                                 ConsensusMap & out
                                ) nogil except +
