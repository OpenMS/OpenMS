from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/PROCESSING/RESAMPLING/LinearResampler.h>" namespace "OpenMS":

    cdef cppclass LinearResampler(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger
        # wrap-doc:
        #  Annotates and filters transitions in a TargetedExperiment
        #  
        #  
        #  :param exp: The input, unfiltered transitions

        LinearResampler() except + nogil 
        LinearResampler(LinearResampler &) except + nogil  # compiler
        void raster(MSSpectrum & input) except + nogil  # wrap-doc:Applies the resampling algorithm to an MSSpectrum
        void rasterExperiment(MSExperiment & input) except + nogil  # wrap-doc:Resamples the data in an MSExperiment
