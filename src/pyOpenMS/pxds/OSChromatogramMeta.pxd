from Types cimport *

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":
    
    cdef cppclass OSChromatogramMeta "OpenSwath::OSChromatogramMeta":

        OSChromatogramMeta() except + nogil  # TODO
        OSChromatogramMeta(OSChromatogramMeta &) except + nogil  # compiler

        size_t index
        libcpp_string id

