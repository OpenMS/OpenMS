from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ScanWindow.h>" namespace "OpenMS":
    
    cdef cppclass ScanWindow "OpenMS::ScanWindow":
        ScanWindow() nogil except +
        ScanWindow(ScanWindow) nogil except +
        DoubleReal begin
        DoubleReal end

