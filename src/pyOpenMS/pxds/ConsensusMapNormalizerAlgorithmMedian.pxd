from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapNormalizerAlgorithmMedian:

        ConsensusMapNormalizerAlgorithmMedian() nogil except +
        ConsensusMapNormalizerAlgorithmMedian(ConsensusMapNormalizerAlgorithmMedian) nogil except + #wrap-ignore

        Size computeMedians(ConsensusMap & input_map, libcpp_vector[double] & medians) nogil except +
        void normalizeMaps(ConsensusMap & input_map, bool shift) nogil except +

