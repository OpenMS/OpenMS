from Types cimport *
from MSExperiment cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerRapid.h>" namespace "OpenMS":
    
    cdef cppclass PeakPickerRapid(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        PeakPickerRapid() nogil except +
        PeakPickerRapid(PeakPickerRapid) nogil except + #wrap-ignore

        void pick(MSSpectrum[Peak1D] & input,
                  MSSpectrum[Peak1D] & output
                 ) nogil except +

        void pickExperiment(MSExperiment[Peak1D, ChromatogramPeak] & input,
                            MSExperiment[Peak1D, ChromatogramPeak] & output
                           ) nogil except +

