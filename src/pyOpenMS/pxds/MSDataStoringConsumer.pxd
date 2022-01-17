from Types cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataStoringConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MSDataStoringConsumer :
        # wrap-doc:
            #   Consumer class that simply stores the data
            #   -----
            #   This class is able to keep spectra and chromatograms passed to it in memory
            #   and the data can be accessed through getData()

        MSDataStoringConsumer() nogil except +
        MSDataStoringConsumer(MSDataStoringConsumer &) nogil except + # compiler

        void setExperimentalSettings(ExperimentalSettings & exp) nogil except + # wrap-doc:Sets experimental settings
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except + # wrap-doc:Sets expected size

        void consumeSpectrum(MSSpectrum & s) nogil except +
        void consumeChromatogram(MSChromatogram & ) nogil except +

        MSExperiment getData() nogil except +
