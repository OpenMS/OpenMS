from Types cimport *
from SignalToNoiseEstimator cimport *
from MSSpectrum cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimatorMeanIterative[Container]:
        # wrap-instances:
        #  SignalToNoiseEstimatorMeanIterative := SignalToNoiseEstimatorMeanIterative[ MSSpectrum ]

        SignalToNoiseEstimatorMeanIterative() except + nogil 
        SignalToNoiseEstimatorMeanIterative(SignalToNoiseEstimatorMeanIterative &) except + nogil  # compiler
        void init(Container & c) except + nogil 
        double getSignalToNoise(Size index) except + nogil 

cdef extern from "<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>" namespace "OpenMS::SignalToNoiseEstimatorMeanIterative<MSSpectrum>":
    
    cdef enum IntensityThresholdCalculation "OpenMS::SignalToNoiseEstimatorMeanIterative<MSSpectrum>::IntensityThresholdCalculation":
        # wrap-attach:
        #    SignalToNoiseEstimatorMeanIterative
        MANUAL
        AUTOMAXBYSTDEV
        AUTOMAXBYPERCENT
