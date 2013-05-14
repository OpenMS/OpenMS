from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>" namespace "OpenMS":

    cdef cppclass PeakPickerCWT(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        PeakPickerCWT()                  nogil except +
        PeakPickerCWT(PeakPickerCWT)   nogil except + #wrap-ignore

        void pick(MSSpectrum[Peak1D] & input, MSSpectrum[Peak1D] & output) nogil except +
        void pickExperiment(MSExperiment[Peak1D, ChromatogramPeak] & input, MSExperiment[Peak1D, ChromatogramPeak] & output) nogil except +
        DoubleReal estimatePeakWidth(MSExperiment[Peak1D, ChromatogramPeak] & input) nogil except +

