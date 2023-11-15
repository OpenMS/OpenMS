# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/Scaler.h>" namespace "OpenMS":

    cdef cppclass Scaler(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        Scaler()       except + nogil 
        Scaler(Scaler &) except + nogil  #wrap-ignore

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 
