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
        #

        ElutionPeakDetection()      nogil except +
        void detectPeaks(MassTrace & in_, libcpp_vector[MassTrace] & out) nogil except +
        void detectPeaks(libcpp_vector[MassTrace] & in_, libcpp_vector[MassTrace] & out) nogil except +
        void filterByPeakWidth(libcpp_vector[MassTrace] & in_, libcpp_vector[MassTrace] & out) nogil except +

        DoubleReal computeMassTraceNoise(MassTrace &) nogil except +
        DoubleReal computeMassTraceSNR(MassTrace &) nogil except +
        DoubleReal computeApexSNR(MassTrace &) nogil except +


