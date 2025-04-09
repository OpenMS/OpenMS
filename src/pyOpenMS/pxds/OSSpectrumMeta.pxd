from Types cimport *

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":
    
    cdef cppclass OSSpectrumMeta "OpenSwath::OSSpectrumMeta":

        OSSpectrumMeta() except + nogil  # TODO
        OSSpectrumMeta(OSSpectrumMeta &) except + nogil  # compiler

        size_t index
        libcpp_string id
        double RT
        int ms_level

