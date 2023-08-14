from Types cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataAggregatingConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MSDataAggregatingConsumer :

        MSDataAggregatingConsumer(MSDataAggregatingConsumer &) except + nogil  # compiler

        void consumeSpectrum(MSSpectrum & s) except + nogil # TODO(whole file)
        void consumeChromatogram(MSChromatogram & ) except + nogil 

        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) except + nogil 
        void setExperimentalSettings(ExperimentalSettings & exp) except + nogil 

