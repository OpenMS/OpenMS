from Types cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimator[Container](DefaultParamHandler,ProgressLogger):
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        SignalToNoiseEstimator() nogil except +
        SignalToNoiseEstimator(SignalToNoiseEstimator) nogil except +
        # void init(Container & c) nogil except +
        # double getSignalToNoise(Size index) nogil except +

