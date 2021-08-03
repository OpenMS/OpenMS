# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>" namespace "OpenMS":

    cdef cppclass BernNorm(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        BernNorm() nogil except +
        BernNorm(BernNorm &) nogil except +

        void filterSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakMap(MSExperiment & exp) nogil except +

