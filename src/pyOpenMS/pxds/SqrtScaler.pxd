# part of the SpectraFilters
from MSSpectrum cimport *
from MSRun cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/PROCESSING/SCALING/SqrtScaler.h>" namespace "OpenMS":

    cdef cppclass SqrtScaler(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        # wrap-doc:
        #  Scales the intensity of peaks to the sqrt

        SqrtScaler() except + nogil  
        SqrtScaler(SqrtScaler &) except + nogil 

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSRun & exp) except + nogil 
