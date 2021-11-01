# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>" namespace "OpenMS":

    cdef cppclass ThresholdMower(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        ThresholdMower() nogil except +
        ThresholdMower(ThresholdMower &) nogil except +

        void filterSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakMap(MSExperiment & exp) nogil except +
