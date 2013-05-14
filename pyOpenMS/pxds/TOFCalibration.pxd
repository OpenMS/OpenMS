from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>" namespace "OpenMS":

    cdef cppclass TOFCalibration(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        TOFCalibration()                  nogil except +
        TOFCalibration(TOFCalibration)   nogil except + #wrap-ignore

        void calibrate(MSExperiment[Peak1D, ChromatogramPeak] & input, MSExperiment[Peak1D, ChromatogramPeak] & output, libcpp_vector[double] & exp_masses) nogil except +
        void pickAndCalibrate(MSExperiment[Peak1D, ChromatogramPeak] & input, MSExperiment[Peak1D, ChromatogramPeak] & output, libcpp_vector[double] & exp_masses) nogil except +
