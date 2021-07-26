from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from InterfaceDataStructures cimport *

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimatorMedianRapid "OpenMS::SignalToNoiseEstimatorMedianRapid":

        SignalToNoiseEstimatorMedianRapid(SignalToNoiseEstimatorMedianRapid &) nogil except + # compiler
        SignalToNoiseEstimatorMedianRapid(double window_length) nogil except +
        NoiseEstimator estimateNoise(shared_ptr[Spectrum]) nogil except +
        NoiseEstimator estimateNoise(shared_ptr[Chromatogram]) nogil except +
        NoiseEstimator estimateNoise(libcpp_vector[ double ] mz_array, libcpp_vector[ double ] int_array) nogil except +

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>" namespace "OpenMS::SignalToNoiseEstimatorMedianRapid":
    
    cdef cppclass NoiseEstimator "OpenMS::SignalToNoiseEstimatorMedianRapid::NoiseEstimator":

        NoiseEstimator() nogil except +
        NoiseEstimator(NoiseEstimator &) nogil except + # compiler
        NoiseEstimator(double nr_windows_, double mz_start_, double win_len_) nogil except +

        int nr_windows
        double mz_start
        double window_length
        libcpp_vector[ double ] result_windows_even
        libcpp_vector[ double ] result_windows_odd

        double get_noise_value(double mz) nogil except +
        double get_noise_even(double mz) nogil except +
        double get_noise_odd(double mz) nogil except +
