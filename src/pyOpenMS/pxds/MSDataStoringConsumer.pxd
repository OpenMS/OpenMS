from Types cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataStoringConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MSDataStoringConsumer :

        MSDataStoringConsumer() nogil except +
        MSDataStoringConsumer(MSDataStoringConsumer) nogil except + #wrap-ignore

        void setExperimentalSettings(ExperimentalSettings & exp) nogil except +
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +

        void consumeSpectrum(MSSpectrum & s) nogil except +
        void consumeChromatogram(MSChromatogram & ) nogil except +

        MSExperiment getData() nogil except +

