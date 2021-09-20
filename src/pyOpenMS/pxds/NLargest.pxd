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
        # wrap-doc:
        #   NLargest removes all but the n largest peaks
        
        NLargest() nogil except + 
        NLargest(NLargest &) nogil except +

        void filterSpectrum(MSSpectrum & spec) nogil except + # wrap-doc:Removes all except spectrums
        void filterPeakSpectrum(MSSpectrum & spec) nogil except + # wrap-doc:Removes all except peak spectrums
        void filterPeakMap(MSExperiment & exp) nogil except + # wrap-doc:Removes all except peak map
