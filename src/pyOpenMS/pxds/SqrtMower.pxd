# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>" namespace "OpenMS":

    cdef cppclass SqrtMower(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        # wrap-doc:
        #  Scales the intensity of peaks to the sqrt

        SqrtMower() except + nogil  
        SqrtMower(SqrtMower &) except + nogil 

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 
