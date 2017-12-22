# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>" namespace "OpenMS":

    cdef cppclass SqrtMower(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        SqrtMower()          nogil except +
        SqrtMower(SqrtMower) nogil except + #wrap-ignore

        void filterSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakMap(MSExperiment & exp) nogil except +

