from Types  cimport *
from smart_ptr cimport shared_ptr
from OpenSwathDataStructures cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>" namespace "OpenSwath":

    cdef cppclass SwathMap:

        SwathMap() nogil except +
        SwathMap(SwathMap) nogil except +
        double lower
        double upper
        bool ms1
        # TODO we would need to support abstract base classes for this ... 
        # OpenSwath::SpectrumAccessPtr sptr;

