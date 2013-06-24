# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>" namespace "OpenMS":

    cdef cppclass WindowMower(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        WindowMower()            nogil except +
        WindowMower(WindowMower) nogil except + #wrap-ignore

        # sliding window version (slower)
        void filterPeakSpectrumForTopNInSlidingWindow(MSSpectrum[Peak1D] & spectrum) nogil except +
        # jumping window version (faster)
        void filterPeakSpectrumForTopNInJumpingWindow(MSSpectrum[Peak1D] & spectrum) nogil except +

        void filterPeakSpectrum(MSSpectrum[Peak1D] & spec) nogil except +
        void filterPeakMap(MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +

