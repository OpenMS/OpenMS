from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/PROCESSING/SMOOTHING/LowessSmoothing.h>" namespace "OpenMS":

    cdef cppclass LowessSmoothing(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        LowessSmoothing() except + nogil 

        void smoothData(libcpp_vector[double] x,
                        libcpp_vector[double] y,
                        libcpp_vector[double] & y_smoothed)      except + nogil  # wrap-doc:Smoothing method that receives x and y coordinates (e.g., RT and intensities) and computes smoothed intensities

