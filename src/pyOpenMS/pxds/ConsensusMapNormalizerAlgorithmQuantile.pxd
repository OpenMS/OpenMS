from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapNormalizerAlgorithmQuantile:

        ConsensusMapNormalizerAlgorithmQuantile() except + nogil 
        # private
        ConsensusMapNormalizerAlgorithmQuantile(ConsensusMapNormalizerAlgorithmQuantile) except + nogil  #wrap-ignore

        void normalizeMaps(ConsensusMap & input_map) except + nogil 

        void resample(libcpp_vector[ double ] & data_in, libcpp_vector[ double ] & data_out, UInt n_resampling_points) except + nogil  # wrap-doc:Resamples data_in and writes the results to data_out
        void extractIntensityVectors(ConsensusMap & map_, libcpp_vector[ libcpp_vector[ double ] ] & out_intensities) except + nogil  # wrap-doc:Extracts the intensities of the features of the different maps
        void setNormalizedIntensityValues(libcpp_vector[ libcpp_vector[ double ] ] & feature_ints, ConsensusMap & map_) except + nogil  # wrap-doc:Writes the intensity values in feature_ints to the corresponding features in map

