# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/NLargest.h>" namespace "OpenMS":

    cdef cppclass NLargest(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        NLargest()            nogil except +
        NLargest(NLargest) nogil except + #wrap-ignore

        void filterSpectrum(MSSpectrum[Peak1D] & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum[Peak1D] & spec) nogil except +
        void filterPeakMap(MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +

