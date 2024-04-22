from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from MassTrace cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

from Types cimport *

from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FEATUREFINDER/ElutionPeakDetection.h>" namespace "OpenMS":

    cdef cppclass ElutionPeakDetection(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #   ProgressLogger
        #   DefaultParamHandler

        ElutionPeakDetection() except + nogil 
        ElutionPeakDetection(ElutionPeakDetection &) except + nogil  # compiler

        void detectPeaks(Kernel_MassTrace & in_,
                         libcpp_vector[Kernel_MassTrace] & out
                         ) except + nogil  # TODO

        void detectPeaks(libcpp_vector[Kernel_MassTrace] & in_,
                         libcpp_vector[Kernel_MassTrace] & out
                        ) except + nogil  # TODO

        void filterByPeakWidth(libcpp_vector[Kernel_MassTrace] & in_,
                               libcpp_vector[Kernel_MassTrace] & out
                              ) except + nogil  # TODO

        double computeMassTraceNoise(Kernel_MassTrace &) except + nogil  # wrap-doc:Compute noise level (as RMSE of the actual signal and the smoothed signal)
        double computeMassTraceSNR(Kernel_MassTrace &) except + nogil  # wrap-doc:Compute the signal to noise ratio (estimated by computeMassTraceNoise)
        double computeApexSNR(Kernel_MassTrace &) except + nogil  # wrap-doc:Compute the signal to noise ratio at the apex (estimated by computeMassTraceNoise)

        void findLocalExtrema(Kernel_MassTrace & , Size & , libcpp_vector[ size_t ] & , libcpp_vector[ size_t ] & ) except + nogil  # TODO

        void smoothData(Kernel_MassTrace & mt, int win_size) except + nogil  # TODO

