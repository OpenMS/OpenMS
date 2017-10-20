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

        void pick(MSSpectrum & input, MSSpectrum & output) nogil except +
        void pickExperiment(MSExperiment & input, MSExperiment & output) nogil except +
        double estimatePeakWidth(MSExperiment & input) nogil except +

