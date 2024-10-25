# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/PROCESSING/FILTERING/NLargest.h>" namespace "OpenMS":

    cdef cppclass NLargest(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        # wrap-doc:
        #  NLargest removes all but the n largest peaks
        
        NLargest() except + nogil  
        NLargest(NLargest &) except + nogil 

        void filterSpectrum(MSSpectrum & spec) except + nogil  # wrap-doc:Keep only n-largest peaks in spectrum
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil  # wrap-doc:Keep only n-largest peaks in spectrum
        void filterPeakMap(MSExperiment & exp) except + nogil  # wrap-doc:Keep only n-largest peaks in each spectrum of a peak map
