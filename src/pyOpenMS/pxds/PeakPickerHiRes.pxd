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

        void pick(MSSpectrum & input,
                  MSSpectrum & output
                 ) nogil except +

        void pickExperiment(MSExperiment & input,
                            MSExperiment & output
                           ) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>" namespace "OpenMS::PeakPickerHiRes":
    
    cdef cppclass PeakBoundary "OpenMS::PeakPickerHiRes::PeakBoundary":
        PeakBoundary() nogil except +
        PeakBoundary(PeakBoundary) nogil except + #wrap-ignore
        double mz_min
        double mz_max
