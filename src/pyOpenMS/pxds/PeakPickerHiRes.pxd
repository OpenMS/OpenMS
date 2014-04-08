from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>" namespace "OpenMS":

    cdef cppclass PeakPickerHiRes(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        PeakPickerHiRes()                  nogil except +
        PeakPickerHiRes(PeakPickerHiRes)   nogil except + #wrap-ignore

        void pick(MSSpectrum[Peak1D] & input,
                  MSSpectrum[Peak1D] & output
                 ) nogil except +

        void pickExperiment(MSExperiment[Peak1D, ChromatogramPeak] & input,
                            MSExperiment[Peak1D, ChromatogramPeak] & output
                           ) nogil except +

