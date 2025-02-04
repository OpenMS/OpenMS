from Types cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataStoringConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MSDataStoringConsumer :
        # wrap-doc:
            #  Consumer class that simply stores the data
            #  
            #  This class is able to keep spectra and chromatograms passed to it in memory
            #  and the data can be accessed through getData()

        MSDataStoringConsumer() except + nogil 
        MSDataStoringConsumer(MSDataStoringConsumer &) except + nogil  # compiler

        void setExperimentalSettings(ExperimentalSettings & exp) except + nogil  # wrap-doc:Sets experimental settings
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) except + nogil  # wrap-doc:Sets expected size

        void consumeSpectrum(MSSpectrum & s) except + nogil 
        void consumeChromatogram(MSChromatogram & ) except + nogil 

        MSExperiment getData() except + nogil 
