from Types cimport *

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":
    
    cdef cppclass OSChromatogramMeta "OpenSwath::OSChromatogramMeta":

        OSChromatogramMeta() nogil except + # TODO
        OSChromatogramMeta(OSChromatogramMeta &) nogil except + # compiler

        size_t index
        libcpp_string id

