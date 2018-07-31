from Types cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>" namespace "OpenMS":
    
    cdef cppclass OSW_ChromExtractParams "OpenMS::ChromExtractParams":

        OSW_ChromExtractParams(OSW_ChromExtractParams) nogil except + #wrap-ignore

        double min_upper_edge_dist
        double mz_extraction_window
        bool ppm
        libcpp_string extraction_function
        double rt_extraction_window
        double extra_rt_extract

