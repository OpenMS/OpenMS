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

        TOFCalibration() nogil except +
        TOFCalibration(TOFCalibration &) nogil except + # compiler

        void calibrate(MSExperiment & input, MSExperiment & output, libcpp_vector[double] & exp_masses) nogil except +
        void pickAndCalibrate(MSExperiment & input, MSExperiment & output, libcpp_vector[double] & exp_masses) nogil except +

        libcpp_vector[ double ]  getML1s() nogil except + # wrap-doc:Returns the first calibration constant
        void setML1s(libcpp_vector[ double ] & ml1s) nogil except +
        libcpp_vector[ double ]  getML2s() nogil except +
        void setML2s(libcpp_vector[ double ] & ml2s) nogil except + # wrap-doc:Returns the second calibration constant
        libcpp_vector[ double ]  getML3s() nogil except +
        void setML3s(libcpp_vector[ double ] & ml3s) nogil except + # wrap-doc:Returns the third calibration constant
