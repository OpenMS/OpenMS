from Types cimport *
from libcpp.set cimport set as libcpp_set
from MSExperiment cimport *
from PeakPickerHiRes cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>" namespace "OpenMS":
    
    cdef cppclass PeakWidthEstimator "OpenMS::PeakWidthEstimator":
        # wrap-doc:
            #   Rough estimation of the peak width at m/z
            #   -----
            #   Based on the peaks of the dataset (peak position & width) and the peak
            #   boundaries as reported by the PeakPickerHiRes, the typical peak width is
            #   estimated for arbitrary m/z using a spline interpolationThis struct can be used to store both peak or feature indices`

        # private
        # PeakWidthEstimator() nogil except +
        PeakWidthEstimator(PeakWidthEstimator &) nogil except + # compiler
        PeakWidthEstimator(MSExperiment exp_picked,
                           libcpp_vector[libcpp_vector[PeakBoundary] ] & boundaries) nogil except +

        double getPeakWidth(double mz) nogil except + # wrap-doc:Returns the estimated peak width at m/z

