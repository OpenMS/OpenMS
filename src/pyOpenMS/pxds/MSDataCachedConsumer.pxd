from MSSpectrum cimport *
from MSChromatogram cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>" namespace "OpenMS":

    cdef cppclass MSDataCachedConsumer:

        MSDataCachedConsumer(String filename) nogil except +
        MSDataCachedConsumer(String filename, bool clear) nogil except +
        MSDataCachedConsumer(MSDataCachedConsumer) nogil except + #wrap-ignore

        void consumeSpectrum(MSSpectrum & s) nogil except +
        void consumeChromatogram(MSChromatogram & c) nogil except +

        void setExperimentalSettings(ExperimentalSettings& exp) nogil except +
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +

