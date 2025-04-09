# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/PROCESSING/FILTERING/WindowMower.h>" namespace "OpenMS":

    cdef cppclass WindowMower(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        WindowMower() except + nogil 
        WindowMower(WindowMower &) except + nogil 

        void filterPeakSpectrumForTopNInSlidingWindow(MSSpectrum & spectrum) except + nogil  # wrap-doc:Sliding window version (slower)
        void filterPeakSpectrumForTopNInJumpingWindow(MSSpectrum & spectrum) except + nogil  # wrap-doc:Jumping window version (faster)

        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 
