from Types cimport *
from SignalToNoiseEstimator cimport *
from MSSpectrum cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimatorMeanIterative[Container]:
        # wrap-instances:
        #   SignalToNoiseEstimatorMeanIterative := SignalToNoiseEstimatorMeanIterative[ MSSpectrum ]

        SignalToNoiseEstimatorMeanIterative() nogil except +
        SignalToNoiseEstimatorMeanIterative(SignalToNoiseEstimatorMeanIterative) nogil except +
        void init(Container & c) nogil except +
        double getSignalToNoise(Size index) nogil except +

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>" namespace "OpenMS::SignalToNoiseEstimatorMeanIterative<MSSpectrum>":
    
    cdef enum IntensityThresholdCalculation "OpenMS::SignalToNoiseEstimatorMeanIterative<MSSpectrum>::IntensityThresholdCalculation":
        # wrap-attach:
        #     SignalToNoiseEstimatorMeanIterative
        MANUAL
        AUTOMAXBYSTDEV
        AUTOMAXBYPERCENT
