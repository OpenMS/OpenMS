# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>" namespace "OpenMS":

    cdef cppclass SpectraMerger(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        SpectraMerger()            nogil except +
        SpectraMerger(SpectraMerger) nogil except + #wrap-ignore

        void mergeSpectraBlockWise(MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +
        void mergeSpectraPrecursors(MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +

