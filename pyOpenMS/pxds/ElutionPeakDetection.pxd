from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from MassTrace cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

from Types cimport *

from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>" namespace "OpenMS":

    cdef cppclass ElutionPeakDetection(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #    ProgressLogger
        #    DefaultParamHandler

        ElutionPeakDetection()      nogil except +

        void detectPeaks(Kernel_MassTrace & in_,
                         libcpp_vector[Kernel_MassTrace] & out
                         ) nogil except +

        void detectPeaks(libcpp_vector[Kernel_MassTrace] & in_,
                         libcpp_vector[Kernel_MassTrace] & out
                        ) nogil except +

        void filterByPeakWidth(libcpp_vector[Kernel_MassTrace] & in_,
                               libcpp_vector[Kernel_MassTrace] & out
                              ) nogil except +

        DoubleReal computeMassTraceNoise(Kernel_MassTrace &) nogil except +
        DoubleReal computeMassTraceSNR(Kernel_MassTrace &) nogil except +
        DoubleReal computeApexSNR(Kernel_MassTrace &) nogil except +

        void findLocalExtrema(Kernel_MassTrace & , Size & , libcpp_vector[ size_t ] & , libcpp_vector[ size_t ] & )

