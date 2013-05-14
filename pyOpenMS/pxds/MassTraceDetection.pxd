from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from MassTrace cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>" namespace "OpenMS":

    cdef cppclass MassTraceDetection(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #    ProgressLogger
        #    DefaultParamHandler

        MassTraceDetection()      nogil except +

        void run(MSExperiment[Peak1D, ChromatogramPeak] & input_map,
                libcpp_vector[MassTrace] & traces
                ) nogil except +

