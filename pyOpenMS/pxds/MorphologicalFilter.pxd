from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>" namespace "OpenMS":

    cdef cppclass MorphologicalFilter(DefaultParamHandler,ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MorphologicalFilter()      nogil except +
        # MorphologicalFilter(MorphologicalFilter)      nogil except + #private
        void filter(MSSpectrum[Peak1D] & spectrum)      nogil except +
        void filterExperiment(MSExperiment[Peak1D,ChromatogramPeak] & exp)      nogil except +

