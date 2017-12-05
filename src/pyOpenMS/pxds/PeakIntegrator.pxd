from DefaultParamHandler cimport *
from String cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>" namespace "OpenMS":

    cdef cppclass PeakIntegrator(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        PeakIntegrator() nogil except +
        PeakIntegrator(PeakIntegrator) nogil except +

        void getDefaultParameters(Param) nogil except +

        void estimateBackground(MSChromatogram, double, double) nogil except +

        void integratePeak(MSChromatogram, double, double) nogil except +

        void calculatePeakShapeMetrics(MSChromatogram, double, double, PeakShapeMetrics_) nogil except +

        double getPeakArea() nogil except +
        double getPeakHeight() nogil except +
        double getPeakApexRT() nogil except +
        double getBackgroundHeight() nogil except +
        double getBackgroundArea() nogil except +