from Types  cimport *
from smart_ptr cimport shared_ptr
from OpenSwathDataStructures cimport *
from ISpectrumAccess cimport *
from smart_ptr cimport shared_ptr

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>" namespace "OpenSwath":

    cdef cppclass SwathMap:

        SwathMap() nogil except +
        SwathMap(SwathMap) nogil except +

        double lower
        double upper
        double center
        bool ms1

        # TODO we would need to support abstract base classes for this ... 
        # OpenSwath::SpectrumAccessPtr sptr;
        ## shared_ptr[ISpectrumAccess] sptr # wrap-ignore

