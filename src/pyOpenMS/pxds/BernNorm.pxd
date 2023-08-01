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
        #   DefaultParamHandler

        BernNorm() except + nogil 
        BernNorm(BernNorm &) except + nogil 

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 

