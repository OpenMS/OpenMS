from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from InterfaceDataStructures cimport *

cdef extern from "<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimatorMedianRapid "OpenMS::SignalToNoiseEstimatorMedianRapid":

        SignalToNoiseEstimatorMedianRapid(SignalToNoiseEstimatorMedianRapid &) except + nogil  # compiler
        SignalToNoiseEstimatorMedianRapid(double window_length) except + nogil 
        NoiseEstimator estimateNoise(shared_ptr[Spectrum]) except + nogil 
        NoiseEstimator estimateNoise(shared_ptr[Chromatogram]) except + nogil 
        NoiseEstimator estimateNoise(libcpp_vector[ double ] mz_array, libcpp_vector[ double ] int_array) except + nogil 

cdef extern from "<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>" namespace "OpenMS::SignalToNoiseEstimatorMedianRapid":
    
    cdef cppclass NoiseEstimator "OpenMS::SignalToNoiseEstimatorMedianRapid::NoiseEstimator":

        NoiseEstimator() except + nogil 
        NoiseEstimator(NoiseEstimator &) except + nogil  # compiler
        NoiseEstimator(double nr_windows_, double mz_start_, double win_len_) except + nogil 

        int nr_windows
        double mz_start
        double window_length
        libcpp_vector[ double ] result_windows_even
        libcpp_vector[ double ] result_windows_odd

        double get_noise_value(double mz) except + nogil 
        double get_noise_even(double mz) except + nogil 
        double get_noise_odd(double mz) except + nogil 
