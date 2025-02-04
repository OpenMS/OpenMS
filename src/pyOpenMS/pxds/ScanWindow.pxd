from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ScanWindow.h>" namespace "OpenMS":
    
    cdef cppclass ScanWindow(MetaInfoInterface) :
        # wrap-inherits:
        #  MetaInfoInterface

        ScanWindow() except + nogil 
        ScanWindow(ScanWindow &) except + nogil 

        double begin
        double end
