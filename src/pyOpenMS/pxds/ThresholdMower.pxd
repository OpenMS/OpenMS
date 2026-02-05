# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/PROCESSING/FILTERING/ThresholdMower.h>" namespace "OpenMS":

    cdef cppclass ThresholdMower(DefaultParamHandler):
        # wrap-doc:
        #  Removes all peaks below an intensity threshold
        # wrap-inherits:
        #   DefaultParamHandler

        ThresholdMower() except + nogil 
        ThresholdMower(ThresholdMower &) except + nogil 

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 
