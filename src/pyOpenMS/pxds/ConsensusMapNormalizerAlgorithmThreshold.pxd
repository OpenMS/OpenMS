from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapNormalizerAlgorithmThreshold:

        ConsensusMapNormalizerAlgorithmThreshold() except + nogil 
        # private
        ConsensusMapNormalizerAlgorithmThreshold(ConsensusMapNormalizerAlgorithmThreshold) except + nogil  # wrap-ignore

        libcpp_vector[double] computeCorrelation(ConsensusMap & input_map,
                                                 double ratio_threshold,
                                                 const String & acc_filter,
                                                 const String & desc_filter) except + nogil  # wrap-doc:Determines the ratio of all maps to the map with the most features

        void normalizeMaps(ConsensusMap & input_map,
                           libcpp_vector[double] & ratios) except + nogil  # wrap-doc:Applies the given ratio to the maps of the consensusMap

