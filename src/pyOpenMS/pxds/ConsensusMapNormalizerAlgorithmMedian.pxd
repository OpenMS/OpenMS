from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *
from ConsensusMapNormalizerAlgorithmThreshold cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>" namespace "OpenMS::ConsensusMapNormalizerAlgorithmMedian":

    cdef enum NormalizationMethod:
        NM_SCALE,
        NM_SHIFT

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>" namespace "OpenMS":
    cdef cppclass ConsensusMapNormalizerAlgorithmMedian:

        ConsensusMapNormalizerAlgorithmMedian() nogil except +
        # private
        ConsensusMapNormalizerAlgorithmMedian(ConsensusMapNormalizerAlgorithmMedian) nogil except + # wrap-ignore

        Size computeMedians(ConsensusMap & input_map,
                            libcpp_vector[double] & medians,
                            const String & acc_filter,
                            const String & desc_filter) nogil except + # wrap-doc:Computes medians of all maps and returns index of map with most features

        void normalizeMaps(ConsensusMap & input_map,
                           NormalizationMethod method,
                           const String & acc_filter,
                           const String & desc_filter) nogil except + # wrap-doc:Normalizes the maps of the consensusMap
