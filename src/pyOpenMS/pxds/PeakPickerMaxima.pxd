from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerMaxima.h>" namespace "OpenMS":

    cdef cppclass PeakPickerMaxima:
        # wrap-doc:
            #   This class implements a fast peak-picking algorithm best suited for
            #   high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
            #   signals of ions with similar mass-to-charge ratios (m/z) exhibit little or
            #   no overlapping and therefore allow for a clear separation. Furthermore, ion
            #   signals tend to show well-defined peak shapes with narrow peak width
            #   -----
            #   This peak-picking algorithm detects ion signals in raw data and
            #   reconstructs the corresponding peak shape by cubic spline interpolation.
            #   Signal detection depends on the signal-to-noise ratio which is adjustable
            #   by the user (see parameter signal_to_noise). A picked peak's m/z and
            #   intensity value is given by the maximum of the underlying peak spline
            #   -----
            #   So far, this peak picker was mainly tested on high resolution data. With
            #   appropriate preprocessing steps (e.g. noise reduction and baseline
            #   subtraction), it might be also applied to low resolution data

        PeakPickerMaxima(double signal_to_noise, double spacing_difference, double sn_window_length) nogil except +
        PeakPickerMaxima(PeakPickerMaxima &) nogil except + # compiler 

        void findMaxima(libcpp_vector[double] mz_array, libcpp_vector[double]
                        int_array, libcpp_vector[PeakCandidate]& pc) nogil except +
            # wrap-doc:
                #   Will find local maxima in raw data
                #   -----
                #   :param mz_array: The array containing m/z values
                #   :param int_array: The array containing intensity values
                #   :param pc: The resulting array containing the peak candidates
                #   :param check_spacings: Check spacing constraints (recommended settings: yes for spectra, no for chromatograms)

        void pick(libcpp_vector[double] mz_array, libcpp_vector[double]
                        int_array, libcpp_vector[PeakCandidate]& pc) nogil except +
          
cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerMaxima.h>" namespace "OpenMS::PeakPickerMaxima":
    
    cdef cppclass PeakCandidate "OpenMS::PeakPickerMaxima::PeakCandidate":

        PeakCandidate() nogil except +
        PeakCandidate(PeakCandidate &) nogil except + # compiler
        int pos
        int left_boundary
        int right_boundary
        double mz_max
        double int_max
