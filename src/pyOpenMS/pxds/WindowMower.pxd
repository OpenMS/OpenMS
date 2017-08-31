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

        void filterPeakSpectrumForTopNInSlidingWindow(MSSpectrum & spectrum) nogil except + # wrap-doc:sliding window version (slower)
        void filterPeakSpectrumForTopNInJumpingWindow(MSSpectrum & spectrum) nogil except + # wrap-doc:jumping window version (faster)

        void filterPeakSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakMap(MSExperiment & exp) nogil except +

