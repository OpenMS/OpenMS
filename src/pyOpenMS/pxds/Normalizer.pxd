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

        Normalizer()           nogil except +
        Normalizer(Normalizer) nogil except + #wrap-ignore

        void filterSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakMap(MSExperiment & exp) nogil except +

