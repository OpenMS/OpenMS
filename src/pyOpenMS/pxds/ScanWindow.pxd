from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ScanWindow.h>" namespace "OpenMS":
    
    cdef cppclass ScanWindow(MetaInfoInterface) :
        # wrap-inherits:
        #  MetaInfoInterface

        ScanWindow() nogil except +
        ScanWindow(ScanWindow) nogil except +

        double begin
        double end

