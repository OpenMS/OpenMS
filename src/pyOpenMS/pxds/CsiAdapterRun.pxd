from Types cimport *
from CsiAdapterIdentification cimport CsiAdapterIdentification

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":
    
    cdef cppclass CsiAdapterRun "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterRun":
        CsiAdapterRun() nogil except +
        CsiAdapterRun(CsiAdapterRun) nogil except + #wrap-ignore
        libcpp_vector[CsiAdapterIdentification] identifications
