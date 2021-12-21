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

        GaussFilter() nogil except + # wrap-doc:This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform profile data
        GaussFilter(GaussFilter &) nogil except + # compiler

        void filter(MSSpectrum & spectrum) nogil except + # wrap-doc:Smoothes an MSSpectrum containing profile data
        void filter(MSChromatogram & chromatogram) nogil except +
        void filterExperiment(MSExperiment & exp) nogil except + # wrap-doc:Smoothes an MSExperiment containing profile data

