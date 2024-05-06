from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>" namespace "OpenMS":

    cdef cppclass SavitzkyGolayFilter(DefaultParamHandler,ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger

        SavitzkyGolayFilter() except + nogil 
        SavitzkyGolayFilter(SavitzkyGolayFilter &) except + nogil  # compiler

        void filter(MSSpectrum & spectrum) except + nogil  # wrap-doc:Removed the noise from an MSSpectrum containing profile data
        void filterExperiment(MSExperiment & exp) except + nogil  # wrap-doc:Removed the noise from an MSExperiment containing profile data
