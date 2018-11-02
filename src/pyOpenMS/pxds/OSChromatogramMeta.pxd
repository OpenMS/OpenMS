from Types cimport *

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":
    
    cdef cppclass OSChromatogramMeta "OpenSwath::OSChromatogramMeta":

        OSChromatogramMeta() nogil except +
        OSChromatogramMeta(OSChromatogramMeta) nogil except + #wrap-ignore

        size_t index
        libcpp_string id

