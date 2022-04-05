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
        PeakPickerSH(PeakPickerSH &) nogil except + # compiler
        void pick(MSSpectrum & input_, MSSpectrum & output, float fWindowWidth) nogil except + # wrap-doc:Applies the peak-picking algorithm to one spectrum
        void pickExperiment(MSExperiment & input_, MSExperiment & output) nogil except + # wrap-doc:Applies the peak-picking algorithm to a map (MSExperiment)

