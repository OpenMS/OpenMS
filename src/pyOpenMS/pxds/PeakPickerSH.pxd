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
        void pick(MSSpectrum & input_, MSSpectrum & output, float fWindowWidth) nogil except +
        void pickExperiment(MSExperiment & input_, MSExperiment & output) nogil except +

