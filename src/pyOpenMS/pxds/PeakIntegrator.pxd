from DefaultParamHandler cimport *
from String cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from Types cimport *
from libcpp.map cimport map as libcpp_map
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>" namespace "OpenMS":

    cdef cppclass PeakIntegrator(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        PeakIntegrator() nogil except +
        PeakIntegrator(PeakIntegrator) nogil except +

        void getDefaultParameters(Param) nogil except +

        void estimateBackground(MSChromatogram, double, double) nogil except +
        void estimateBackground(MSSpectrum, double, double) nogil except +

        void integratePeak(MSChromatogram, double, double) nogil except +
        void integratePeak(MSSpectrum, double, double) nogil except +

        void calculatePeakShapeMetrics(MSChromatogram, double, double) nogil except +
        void calculatePeakShapeMetrics(MSSpectrum, double, double) nogil except +

        double getPeakArea() nogil except +
        double getPeakHeight() nogil except +
        double getPeakApexRT() nogil except +
        double getBackgroundHeight() nogil except +
        double getBackgroundArea() nogil except +
        double getWidthAt5() nogil except +
        double getWidthAt10() nogil except +
        double getWidthAt50() nogil except +
        double getStartTimeAt5() nogil except +
        double getStartTimeAt10() nogil except +
        double getStartTimeAt50() nogil except +
        double getEndTimeAt5() nogil except +
        double getEndTimeAt10() nogil except +
        double getEndTimeAt50() nogil except +
        double getTotalWidth() nogil except +
        double getTailingFactor() nogil except +
        double getAsymmetryFactor() nogil except +
        double getBaselineDeltaToHeight() nogil except +
        double getSlopeOfBaseline() nogil except +
        Int getPointsAcrossBaseline() nogil except +
        Int getPointsAcrossHalfHeight() nogil except +

        libcpp_map[ String, double ] getPeakShapeMetrics() nogil except +
