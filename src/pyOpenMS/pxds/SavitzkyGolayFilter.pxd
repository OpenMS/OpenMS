from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>" namespace "OpenMS":

    cdef cppclass SavitzkyGolayFilter(DefaultParamHandler,ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        SavitzkyGolayFilter()      nogil except +
        SavitzkyGolayFilter(SavitzkyGolayFilter)      nogil except +

        void filter(MSSpectrum & spectrum)      nogil except +
        void filterExperiment(MSExperiment & exp)      nogil except +


