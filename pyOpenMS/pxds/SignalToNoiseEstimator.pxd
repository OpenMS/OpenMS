from Types cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>" namespace "OpenMS":
    
    cdef cppclass SignalToNoiseEstimator[Container](DefaultParamHandler,ProgressLogger):
        # wrap-ignore
        # ABSTRACT class
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        SignalToNoiseEstimator() nogil except +
        SignalToNoiseEstimator(SignalToNoiseEstimator) nogil except +
        # void init(PeakIterator & it_begin, PeakIterator & it_end) nogil except +
        # void init(Container & c) nogil except +
        # double getSignalToNoise(PeakIterator & data_point) nogil except +
        # double getSignalToNoise(Peak1D & data_point) nogil except +

