from Types cimport *
from SiriusAdapterIdentification cimport SiriusAdapterIdentification

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterRun "OpenMS::SiriusMzTabWriter::SiriusAdapterRun":
        SiriusAdapterRun() nogil except + 
        SiriusAdapterRun(SiriusAdapterRun) nogil except + #wrap-ignore
        libcpp_vector[SiriusAdapterIdentification] identifications
