from ConvexHull2D cimport *
from DefaultParamHandler cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>" namespace "OpenMS":

    cdef cppclass PeakIntegrator(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #   Compute the area, background and shape metrics of a peak
        #   -----
        #   The area computation is performed in integratePeak() and it supports
        #   integration by simple sum of the intensity, integration by Simpson's rule
        #   implementations for an odd number of unequally spaced points or integration
        #   by the trapezoid rule
        #   -----
        #   The background computation is performed in estimateBackground() and it
        #   supports three different approaches to baseline correction, namely
        #   computing a rectangular shape under the peak based on the minimum value of
        #   the peak borders (vertical_division_min), a rectangular shape based on the
        #   maximum value of the beak borders (vertical_division_max) or a trapezoidal
        #   shape based on a straight line between the peak borders (base_to_base)
        #   -----
        #   Peak shape metrics are computed in calculatePeakShapeMetrics() and multiple
        #   metrics are supported
        #   -----
        #   The containers supported by the methods are MSChromatogram and MSSpectrum

        PeakIntegrator() nogil except +
        PeakIntegrator(PeakIntegrator &) nogil except + # compiler

        void getDefaultParameters(Param) nogil except +

        PI_PeakArea integratePeak(MSChromatogram& chromatogram, double left, double right) nogil except +
        PI_PeakArea integratePeak(MSSpectrum& spectrum, double left, double right) nogil except +

        PI_PeakBackground estimateBackground(MSChromatogram& chromatogram, double left, double right, double peak_apex_pos) nogil except +
        PI_PeakBackground estimateBackground(MSSpectrum& spectrum, double left, double right, double peak_apex_pos) nogil except +

        PI_PeakShapeMetrics calculatePeakShapeMetrics(MSChromatogram& chromatogram, double left, double right, double peak_height, double peak_apex_pos) nogil except +
        PI_PeakShapeMetrics calculatePeakShapeMetrics(MSSpectrum& spectrum, double left, double right, double peak_height, double peak_apex_pos) nogil except +


cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>" namespace "OpenMS::PeakIntegrator":

    cdef cppclass PI_PeakArea "OpenMS::PeakIntegrator::PeakArea":

        PI_PeakArea() nogil except +
        PI_PeakArea(PI_PeakArea &) nogil except + # compiler

        double area
        double height
        double apex_pos
        libcpp_vector[DPosition2] hull_points

    cdef cppclass PI_PeakBackground "OpenMS::PeakIntegrator::PeakBackground":

        PI_PeakBackground() nogil except +
        PI_PeakBackground(PI_PeakBackground &) nogil except + # compiler

        double area
        double height

    cdef cppclass PI_PeakShapeMetrics "OpenMS::PeakIntegrator::PeakShapeMetrics":

        PI_PeakShapeMetrics() nogil except +
        PI_PeakShapeMetrics(PI_PeakShapeMetrics &) nogil except + # compiler

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
