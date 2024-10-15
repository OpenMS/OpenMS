from Types cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimator[Container](DefaultParamHandler,ProgressLogger):
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        # wrap-doc:
        #  This class represents the abstract base class of a signal to noise estimator
        
        SignalToNoiseEstimator() except + nogil 
        SignalToNoiseEstimator(SignalToNoiseEstimator &) except + nogil 
        # void init(Container & c) except + nogil 
        # double getSignalToNoise(Size index) except + nogil 
