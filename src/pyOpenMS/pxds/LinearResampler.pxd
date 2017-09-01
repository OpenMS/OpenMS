from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>" namespace "OpenMS":

    cdef cppclass LinearResampler(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        LinearResampler()                  nogil except +
        LinearResampler(LinearResampler)   nogil except + #wrap-ignore
        void raster(MSSpectrum & input) nogil except +
        void rasterExperiment(MSExperiment & input) nogil except +

