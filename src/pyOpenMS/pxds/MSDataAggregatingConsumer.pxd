from Types cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataAggregatingConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MSDataAggregatingConsumer :

        MSDataAggregatingConsumer(MSDataAggregatingConsumer &) nogil except + # compiler

        void consumeSpectrum(MSSpectrum & s) nogil except +# TODO(whole file)
        void consumeChromatogram(MSChromatogram & ) nogil except +

        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +
        void setExperimentalSettings(ExperimentalSettings & exp) nogil except +

