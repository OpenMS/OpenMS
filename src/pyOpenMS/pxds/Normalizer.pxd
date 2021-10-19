# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>" namespace "OpenMS":

    cdef cppclass Normalizer(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler
        # wrap-doc:
        #   Normalizes the peak intensities spectrum-wise

        Normalizer() nogil except +

        Normalizer(Normalizer) nogil except +

        void filterSpectrum(MSSpectrum & spec) nogil except + # wrap-doc:Normalizes the spectrum

        void filterPeakSpectrum(MSSpectrum & spec) nogil except + # wrap-doc:Normalizes the peak spectrum

        void filterPeakMap(MSExperiment & exp) nogil except + # wrap-doc:Normalizes the peak map
