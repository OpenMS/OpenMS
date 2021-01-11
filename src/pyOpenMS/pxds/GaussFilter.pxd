from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/SMOOTHING/GaussFilter.h>" namespace "OpenMS":

    cdef cppclass GaussFilter(DefaultParamHandler,ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        GaussFilter() nogil except +
        GaussFilter(GaussFilter) nogil except +

        void filter(MSSpectrum & spectrum) nogil except +
        void filter(MSChromatogram & chromatogram) nogil except +
        void filterExperiment(MSExperiment & exp) nogil except +

