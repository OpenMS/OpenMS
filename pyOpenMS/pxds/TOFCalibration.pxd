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

        libcpp_vector[ double ]  getML1s()
        void setML1s(libcpp_vector[ double ] & ml1s)
        libcpp_vector[ double ]  getML2s()
        void setML2s(libcpp_vector[ double ] & ml2s)
        libcpp_vector[ double ]  getML3s()
        void setML3s(libcpp_vector[ double ] & ml3s)

