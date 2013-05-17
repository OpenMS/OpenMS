from Types cimport *
from SignalToNoiseEstimator cimport *

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimatorMeanIterative[Container](DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        # wrap-instances:
        #   SignalToNoiseEstimator := SignalToNoiseEstimator[ MSSpectrum[Peak1D] ]
        SignalToNoiseEstimatorMeanIterative() nogil except +
        SignalToNoiseEstimatorMeanIterative(SignalToNoiseEstimatorMeanIterative) nogil except +
        # void init(PeakIterator & it_begin, PeakIterator & it_end) nogil except +
        void init(Container & c) nogil except +
        # double getSignalToNoise(PeakIterator & data_point) nogil except +
        # double getSignalToNoise(Peak1D & data_point) nogil except +

# FAILS TODO
# #cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>" namespace "OpenMS::SignalToNoiseEstimatorMeanIterative":
# #    cdef enum IntensityThresholdCalculation "OpenMS::SignalToNoiseEstimatorMeanIterative::IntensityThresholdCalculation":
# #        #wrap-attach:
# #        #    SignalToNoiseEstimatorMeanIterative
# #        MANUAL
# #        AUTOMAXBYSTDEV
# #        AUTOMAXBYPERCENT
# #
