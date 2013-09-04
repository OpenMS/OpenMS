from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerIterative.h>" namespace "OpenMS":

    cdef cppclass PeakPickerIterative(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        PeakPickerIterative() nogil except +
        PeakPickerIterative(PeakPickerIterative) nogil except + #wrap-ignore

        void pick(MSSpectrum[Peak1D] & input,
                  MSSpectrum[Peak1D] & output
                 ) nogil except +

        void pickExperiment(MSExperiment[Peak1D, ChromatogramPeak] & input,
                            MSExperiment[Peak1D, ChromatogramPeak] & output
                           ) nogil except +

