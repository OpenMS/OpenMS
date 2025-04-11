from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSRun cimport *
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
        void filterExperiment(MSRun & exp) except + nogil  # wrap-doc:Removed the noise from an MSRun containing profile data
