from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from SiriusAdapterHit cimport *
from MzTab cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMzTabWriter "OpenMS::SiriusMzTabWriter":
        SiriusMzTabWriter() nogil except + 
        SiriusMzTabWriter(SiriusMzTabWriter) nogil except + #wrap-ignore
        String extract_scan_index(const String & path) nogil except +
        void read(libcpp_vector[String] sirius_output_paths, String original_input_mzml, Size top_n_hits, MzTab & result) nogil except +


cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterRun "OpenMS::SiriusMzTabWriter::SiriusAdapterRun":
        SiriusAdapterRun() nogil except + 
        SiriusAdapterRun(SiriusAdapterRun) nogil except + #wrap-ignore
        libcpp_vector[ SiriusAdapterIdentification ] identifications


cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterIdentification "OpenMS::SiriusMzTabWriter::SiriusAdapterIdentification":
        SiriusAdapterIdentification() nogil except +
        SiriusAdapterIdentification(SiriusAdapterIdentification) nogil except + # wrap-ignore

        String scan_index
        libcpp_vector[ SiriusAdapterHit ] hits

