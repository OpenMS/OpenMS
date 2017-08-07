from Types cimport *
from String cimport *
from MzTab cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass CsiFingerIdMzTabWriter "OpenMS::CsiFingerIdMzTabWriter":
        CsiFingerIdMzTabWriter() nogil except + 
        CsiFingerIdMzTabWriter(CsiFingerIdMzTabWriter) nogil except + #wrap-ignore
        void read(libcpp_vector[ String ] & paths, Size number, MzTab & result) nogil except +

