from Types cimport *
from SiriusAdapterIdentification cimport SiriusAdapterIdentification

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterRun "OpenMS::SiriusMzTabWriter::SiriusAdapterRun":
        SiriusAdapterRun() except + nogil  
        SiriusAdapterRun(SiriusAdapterRun &) except + nogil  # compiler

        libcpp_vector[SiriusAdapterIdentification] identifications
