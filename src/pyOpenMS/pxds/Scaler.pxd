# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/Scaler.h>" namespace "OpenMS":

    cdef cppclass Scaler(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        Scaler()       nogil except +
        Scaler(Scaler) nogil except + #wrap-ignore

        void filterSpectrum(MSSpectrum[Peak1D] & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum[Peak1D] & spec) nogil except +
        void filterPeakMap(MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +

