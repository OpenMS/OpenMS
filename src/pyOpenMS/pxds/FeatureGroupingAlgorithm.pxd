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
        # Abstract Class
        void transferSubelements(libcpp_vector[ConsensusMap] maps,
                                 ConsensusMap & out
                                ) nogil except + # wrap-doc:Transfers subelements (grouped features) from input consensus maps to the result consensus map

        # since this is a base class, cannot have overloaded methods
        # void group(libcpp_vector[ FeatureMap ] & maps, ConsensusMap & out)
        # void group(libcpp_vector[ ConsensusMap ] & maps, ConsensusMap & out)
        void registerChildren() nogil except + # wrap-doc:Register all derived classes in this method
