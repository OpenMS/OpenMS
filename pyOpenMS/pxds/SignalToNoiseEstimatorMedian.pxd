from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>" namespace "OpenMS":

    cdef cppclass SignalToNoiseEstimatorMedian[SpectrumT]:
        # wrap-instances:
        #   SignalToNoiseEstimatorMedian := SignalToNoiseEstimatorMedian[MSSpectrum[Peak1D]]

        SignalToNoiseEstimatorMedian() nogil except +
        SignalToNoiseEstimatorMedian(SignalToNoiseEstimatorMedian) nogil except +

        void init(MSSpectrum[Peak1D] & spectrum) nogil except +
        void getSignalToNoise(Peak1D & data_point) nogil except +

        # Functions for SignalToNoiseEstimatorMedianChrom[MSChromatogram[ChromatogramPeak]]
        # use wrap-ignore because autowrap cannot handle them at the moment
        # see addons/SignalToNoiseEstimatorMedianChrom.pyx for the implementation
        void init(MSChromatogram[ChromatogramPeak] & spectrum) nogil except + #wrap-ignore
        void getSignalToNoise(ChromatogramPeak & data_point) nogil except + #wrap-ignore

