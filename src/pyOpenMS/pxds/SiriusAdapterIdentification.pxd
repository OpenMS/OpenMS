from Types cimport *
from SiriusAdapterHit cimport *
from String cimport *
from StringList cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":

    cdef cppclass SiriusAdapterIdentification "OpenMS::SiriusMzTabWriter::SiriusAdapterIdentification":
        SiriusAdapterIdentification() nogil except +
        SiriusAdapterIdentification(SiriusAdapterIdentification &) nogil except + # compiler

        double mz
        double rt
        StringList native_ids
        int scan_index
        int scan_number
        String feature_id
        libcpp_vector[ SiriusAdapterHit ] hits
