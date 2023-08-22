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
        #   DefaultParamHandler
        #   ProgressLogger

        TOFCalibration() except + nogil 
        TOFCalibration(TOFCalibration &) except + nogil  # compiler

        void calibrate(MSExperiment & input, MSExperiment & output, libcpp_vector[double] & exp_masses) except + nogil 
        void pickAndCalibrate(MSExperiment & input, MSExperiment & output, libcpp_vector[double] & exp_masses) except + nogil 

        libcpp_vector[ double ]  getML1s() except + nogil  # wrap-doc:Returns the first calibration constant
        void setML1s(libcpp_vector[ double ] & ml1s) except + nogil 
        libcpp_vector[ double ]  getML2s() except + nogil 
        void setML2s(libcpp_vector[ double ] & ml2s) except + nogil  # wrap-doc:Returns the second calibration constant
        libcpp_vector[ double ]  getML3s() except + nogil 
        void setML3s(libcpp_vector[ double ] & ml3s) except + nogil  # wrap-doc:Returns the third calibration constant
