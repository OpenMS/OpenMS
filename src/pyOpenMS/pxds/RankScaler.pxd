# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/PROCESSING/SCALING/RankScaler.h>" namespace "OpenMS":

    cdef cppclass RankScaler(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        RankScaler()       except + nogil 
        RankScaler(RankScaler &) except + nogil  #wrap-ignore

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 
