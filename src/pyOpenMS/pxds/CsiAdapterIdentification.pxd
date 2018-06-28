from Types cimport *
from CsiAdapterHit cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":
    
    cdef cppclass CsiAdapterIdentification "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterIdentification":
        CsiAdapterIdentification() nogil except +
        CsiAdapterIdentification(CsiAdapterIdentification) nogil except + #wrap-ignore
        int scan_index
        int scan_number
        String feature_id
        libcpp_vector[ CsiAdapterHit ] hits

