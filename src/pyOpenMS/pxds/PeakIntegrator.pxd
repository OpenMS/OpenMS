from ConvexHull2D cimport *
from DefaultParamHandler cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>" namespace "OpenMS":

    cdef cppclass PeakIntegrator(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        PeakIntegrator() nogil except +
        PeakIntegrator(PeakIntegrator) nogil except +

        void getDefaultParameters(Param) nogil except +

        PI_PeakArea integratePeak(MSChromatogram chromatogram, double left, double right) nogil except +
        PI_PeakArea integratePeak(MSSpectrum spectrum, double left, double right) nogil except +

        PI_PeakBackground estimateBackground(MSChromatogram chromatogram, double left, double right, double peak_apex_pos) nogil except +
        PI_PeakBackground estimateBackground(MSSpectrum spectrum, double left, double right, double peak_apex_pos) nogil except +

        PI_PeakShapeMetrics calculatePeakShapeMetrics(MSChromatogram chromatogram, double left, double right, double peak_height, double peak_apex_pos) nogil except +
        PI_PeakShapeMetrics calculatePeakShapeMetrics(MSSpectrum spectrum, double left, double right, double peak_height, double peak_apex_pos) nogil except +

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>" namespace "OpenMS::PeakIntegrator":

    cdef cppclass PI_PeakArea "OpenMS::PeakIntegrator::PeakArea":

        PI_PeakArea() nogil except +
        PI_PeakArea(PI_PeakArea) nogil except +

        double area
        double height
        double apex_pos
        libcpp_vector[DPosition2] hull_points

    cdef cppclass PI_PeakBackground "OpenMS::PeakIntegrator::PeakBackground":

        PI_PeakBackground() nogil except +
        PI_PeakBackground(PI_PeakBackground) nogil except +

        double area
        double height

    cdef cppclass PI_PeakShapeMetrics "OpenMS::PeakIntegrator::PeakShapeMetrics":

        PI_PeakShapeMetrics() nogil except +
        PI_PeakShapeMetrics(PI_PeakShapeMetrics) nogil except +

        double width_at_5
        double width_at_10
        double width_at_50
        double start_position_at_5
        double start_position_at_10
        double start_position_at_50
        double end_position_at_5
        double end_position_at_10
        double end_position_at_50
        double total_width
        double tailing_factor
        double asymmetry_factor
        double slope_of_baseline
        double baseline_delta_2_height
        Int points_across_baseline
        Int points_across_half_height
