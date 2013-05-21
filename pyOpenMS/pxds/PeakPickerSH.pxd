from Types cimport *
from MSExperiment cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerSH.h>" namespace "OpenMS":
    
    cdef cppclass PeakPickerSH(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        PeakPickerSH() nogil except +
        PeakPickerSH(PeakPickerSH) nogil except + #wrap-ignore
        # TODO OpenMS API problem: even though a template, the funciton is not in the header! -> it cannot be called 
        # void pick(MSSpectrum[ Peak1D ] & input_, MSSpectrum[ Peak1D ] & output, float fWindowWidth) nogil except +
        void pickExperiment(MSExperiment[Peak1D, ChromatogramPeak] & input_, MSExperiment[Peak1D, ChromatogramPeak] & output) nogil except +

