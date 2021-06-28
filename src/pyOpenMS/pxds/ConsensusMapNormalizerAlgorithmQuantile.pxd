from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapNormalizerAlgorithmQuantile:

        ConsensusMapNormalizerAlgorithmQuantile() nogil except +
        ConsensusMapNormalizerAlgorithmQuantile(ConsensusMapNormalizerAlgorithmQuantile) nogil except + #wrap-ignore

        void normalizeMaps(ConsensusMap & input_map) nogil except +

        void resample(libcpp_vector[ double ] & data_in, libcpp_vector[ double ] & data_out, UInt n_resampling_points) nogil except + # wrap-doc:resamples data_in and writes the results to data_out
        void extractIntensityVectors(ConsensusMap & map_, libcpp_vector[ libcpp_vector[ double ] ] & out_intensities) nogil except + # wrap-doc:extracts the intensities of the features of the different maps
        void setNormalizedIntensityValues(libcpp_vector[ libcpp_vector[ double ] ] & feature_ints, ConsensusMap & map_) nogil except + # wrap-doc:writes the intensity values in feature_ints to the corresponding features in map

