from Types cimport *
from libcpp.set cimport set as libcpp_set
from MSExperiment cimport *
from PeakPickerHiRes cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>" namespace "OpenMS":
    
    cdef cppclass PeakWidthEstimator "OpenMS::PeakWidthEstimator":

        # PeakWidthEstimator() nogil except + # default constructor is private
        PeakWidthEstimator(PeakWidthEstimator) nogil except + #wrap-ignore
        PeakWidthEstimator(MSExperiment exp_picked,
                           libcpp_vector[libcpp_vector[PeakBoundary] ] & boundaries) nogil except +

        double getPeakWidth(double mz) nogil except +

