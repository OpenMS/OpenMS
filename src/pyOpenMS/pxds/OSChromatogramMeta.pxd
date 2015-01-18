from Types cimport *
from OpenSwathAlgoConfig cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":
    
    cdef cppclass OSChromatogramMeta "OpenSwath::OSChromatogramMeta":
        OSChromatogramMeta() nogil except +
        OSChromatogramMeta(OSChromatogramMeta) nogil except + #wrap-ignore
        # NAMESPACE # std::size_t index
        libcpp_string id

