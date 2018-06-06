from Types cimport *
from smart_ptr cimport shared_ptr
from OpenSwathDataStructures cimport *
from ISpectrumAccess cimport *

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>" namespace "OpenSwath":

    cdef cppclass SwathMap:

        SwathMap() nogil except +
        SwathMap(SwathMap) nogil except +

        double lower
        double upper
        double center
        bool ms1

        # COMMENT: access through
        # COMMENT:  - getSpectrumPtr
        # COMMENT:  - setSpectrumPtr
        shared_ptr[ISpectrumAccess] sptr # wrap-ignore

