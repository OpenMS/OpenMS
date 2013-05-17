from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from DefaultParamHandler cimport *
from BaseFeature cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>" namespace "OpenMS":
    
    cdef cppclass FeatureDistance(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        FeatureDistance(FeatureDistance) nogil except + #wrap-ignore
        # TODO  is static const -> no setters please
        # DoubleReal infinity
        FeatureDistance(DoubleReal max_intensity, bool force_constraints) nogil except +
        # libcpp_pair[ bool, DoubleReal ] operator()(BaseFeature & left, BaseFeature & right) nogil except +

