from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from DefaultParamHandler cimport *
from BaseFeature cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>" namespace "OpenMS":
    
    cdef cppclass FeatureDistance(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        FeatureDistance(double max_intensity, bool force_constraints) except + nogil 
        FeatureDistance(FeatureDistance &) except + nogil  # compiler

        # TODO  is static const -> no setters please
        # double infinity
        # libcpp_pair[ bool, double ] operator()(BaseFeature & left, BaseFeature & right) except + nogil 

