from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>" namespace "OpenMS":

    cdef cppclass GaussFilter(DefaultParamHandler,ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger

        GaussFilter() except + nogil  # wrap-doc:This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform profile data
        GaussFilter(GaussFilter &) except + nogil  # compiler

        void filter(MSSpectrum & spectrum) except + nogil  # wrap-doc:Smoothes an MSSpectrum containing profile data
        void filter(MSChromatogram & chromatogram) except + nogil 
        void filterExperiment(MSExperiment & exp) except + nogil  # wrap-doc:Smoothes an MSExperiment containing profile data

