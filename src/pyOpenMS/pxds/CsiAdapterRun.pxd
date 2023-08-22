from Types cimport *
from CsiAdapterIdentification cimport CsiAdapterIdentification

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":
    
    cdef cppclass CsiAdapterRun "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterRun":
        CsiAdapterRun() except + nogil 
        CsiAdapterRun(CsiAdapterRun &) except + nogil  # compiler
        libcpp_vector[CsiAdapterIdentification] identifications
