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
        #   SignalToNoiseEstimatorMedian := SignalToNoiseEstimatorMedian[MSSpectrum]

        SignalToNoiseEstimatorMedian() nogil except +
        SignalToNoiseEstimatorMedian(SignalToNoiseEstimatorMedian) nogil except +

        void init(MSSpectrum & spectrum) nogil except +
        double getSignalToNoise(Size index) nogil except +

        double getSparseWindowPercent() nogil except +
        double getHistogramRightmostPercent() nogil except +

        # Functions for SignalToNoiseEstimatorMedianChrom[MSChromatogram]
        # use wrap-ignore because autowrap cannot handle them at the moment
        # see addons/SignalToNoiseEstimatorMedianChrom.pyx for the implementation
        void init(MSChromatogram & spectrum) nogil except + #wrap-ignore

cdef extern from "<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>" namespace "OpenMS::SignalToNoiseEstimatorMedian":

    cdef enum IntensityThresholdCalculation "OpenMS::SignalToNoiseEstimatorMedianChrom::IntensityThresholdCalculation":
        # wrap-attach:
        #     SignalToNoiseEstimatorMedian
        MANUAL
        AUTOMAXBYSTDEV
        AUTOMAXBYPERCENT
