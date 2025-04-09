from Types cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from MassTrace cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FEATUREFINDER/MassTraceDetection.h>" namespace "OpenMS":

    cdef cppclass MassTraceDetection(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #   ProgressLogger
        #   DefaultParamHandler

        MassTraceDetection() except + nogil 
        MassTraceDetection(MassTraceDetection &) except + nogil  # compiler

        void run(MSExperiment & input_map,
                libcpp_vector[Kernel_MassTrace] & traces,
                Size max_traces) except + nogil 
