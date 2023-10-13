from Types cimport *
from CsiAdapterHit cimport *
from String cimport *
from StringList cimport *
from libcpp.vector cimport vector as libcpp_vector


cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":
    
    cdef cppclass CsiAdapterIdentification "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterIdentification":
        CsiAdapterIdentification() except + nogil 
        CsiAdapterIdentification(CsiAdapterIdentification& ) except + nogil  # compiler

        double mz
        double rt
        StringList native_ids
        int scan_index
        int scan_number
        String feature_id
        libcpp_vector[CsiAdapterHit] hits
