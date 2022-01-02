from Types cimport *

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":
    
    cdef cppclass OSSpectrumMeta "OpenSwath::OSSpectrumMeta":

        OSSpectrumMeta() nogil except + # TODO
        OSSpectrumMeta(OSSpectrumMeta &) nogil except + # compiler

        size_t index
        libcpp_string id
        double RT
        int ms_level

