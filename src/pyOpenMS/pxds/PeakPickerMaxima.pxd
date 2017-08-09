from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerMaxima.h>" namespace "OpenMS":

    cdef cppclass PeakPickerMaxima:

        PeakPickerMaxima() nogil except + # wrap-ignore
        PeakPickerMaxima(PeakPickerMaxima &) nogil except +
        PeakPickerMaxima(double signal_to_noise, double spacing_difference, double sn_window_length) nogil except +

        void findMaxima(libcpp_vector[double] mz_array, libcpp_vector[double]
                        int_array, libcpp_vector[PeakCandidate]& pc) nogil except +

        void pick(libcpp_vector[double] mz_array, libcpp_vector[double]
                        int_array, libcpp_vector[PeakCandidate]& pc) nogil except +
          
cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerMaxima.h>" namespace "OpenMS::PeakPickerMaxima":
    
    cdef cppclass PeakCandidate "OpenMS::PeakPickerMaxima::PeakCandidate":

        PeakCandidate() nogil except +
        PeakCandidate(PeakCandidate &) nogil except +
        int pos
        int left_boundary
        int right_boundary
        double mz_max
        double int_max
