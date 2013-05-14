from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapNormalizerAlgorithmQuantile:

        ConsensusMapNormalizerAlgorithmQuantile() nogil except +
        ConsensusMapNormalizerAlgorithmQuantile(ConsensusMapNormalizerAlgorithmQuantile) nogil except + #wrap-ignore

        void normalizeMaps(ConsensusMap & input_map) nogil except +

