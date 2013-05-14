from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>" namespace "OpenMS":

    cdef cppclass LowessSmoothing(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        LowessSmoothing()      nogil except +
        LowessSmoothing(LowessSmoothing)      nogil except +

        void smoothData(libcpp_vector[double] x,
                        libcpp_vector[double] y,
                        libcpp_vector[double] & y_smoothed)      nogil except +

