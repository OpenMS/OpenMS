from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>" namespace "OpenMS":

    cdef cppclass SignalToNoiseEstimatorMedian[SpectrumT]:
        # wrap-instances:
        #  SignalToNoiseEstimatorMedian := SignalToNoiseEstimatorMedian[MSSpectrum]

        SignalToNoiseEstimatorMedian() except + nogil 
        SignalToNoiseEstimatorMedian(SignalToNoiseEstimatorMedian &) except + nogil  # compiler

        void init(MSSpectrum & spectrum) except + nogil 
        double getSignalToNoise(Size index) except + nogil 

        double getSparseWindowPercent() except + nogil 
        double getHistogramRightmostPercent() except + nogil 

        # Functions for SignalToNoiseEstimatorMedianChrom[MSChromatogram]
        # use wrap-ignore because autowrap cannot handle them at the moment
        # see addons/SignalToNoiseEstimatorMedianChrom.pyx for the implementation
        void init(MSChromatogram & spectrum) except + nogil  #wrap-ignore

cdef extern from "<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>" namespace "OpenMS::SignalToNoiseEstimatorMedian":

    cdef enum IntensityThresholdCalculation "OpenMS::SignalToNoiseEstimatorMedianChrom::IntensityThresholdCalculation":
        # wrap-attach:
        #    SignalToNoiseEstimatorMedian
        MANUAL
        AUTOMAXBYSTDEV
        AUTOMAXBYPERCENT
