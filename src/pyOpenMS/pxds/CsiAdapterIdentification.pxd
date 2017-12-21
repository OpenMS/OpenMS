from Types cimport *
from CsiAdapterHit cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":
    
    cdef cppclass CsiAdapterIdentification "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterIdentification":
        CsiAdapterIdentification() nogil except +
        CsiAdapterIdentification(CsiAdapterIdentification) nogil except + #wrap-ignore
        String scan_index
        libcpp_vector[ CsiAdapterHit ] hits

