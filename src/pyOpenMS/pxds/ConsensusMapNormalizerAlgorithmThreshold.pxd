from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapNormalizerAlgorithmThreshold:

        ConsensusMapNormalizerAlgorithmThreshold() nogil except +
        ConsensusMapNormalizerAlgorithmThreshold(ConsensusMapNormalizerAlgorithmThreshold) nogil except + #wrap-ignore

        libcpp_vector[double] computeCorrelation(ConsensusMap & input_map,
                                                 double ratio_threshold,
                                                 const String & acc_filter,
                                                 const String & desc_filter) nogil except +

        void normalizeMaps(ConsensusMap & input_map,
                           libcpp_vector[double] & ratios) nogil except +

