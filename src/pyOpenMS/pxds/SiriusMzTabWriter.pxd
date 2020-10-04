from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from SiriusAdapterHit cimport *
from MzTab cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMzTabWriter "OpenMS::SiriusMzTabWriter":
        SiriusMzTabWriter() nogil except + 
        SiriusMzTabWriter(SiriusMzTabWriter) nogil except + #wrap-ignore
        int extractScanIndex(const String& path) nogil except +
        int extractScanNumber(const String& path) nogil except +
        String extractFeatureId(const String& path) nogil except +
        void read(libcpp_vector[String] sirius_output_paths, String original_input_mzml, Size top_n_hits, MzTab & result) nogil except +

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterRun "OpenMS::SiriusMzTabWriter::SiriusAdapterRun":
        SiriusAdapterRun() nogil except + 
        SiriusAdapterRun(SiriusAdapterRun) nogil except + #wrap-ignore
        libcpp_vector[SiriusAdapterIdentification] identifications


cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterIdentification "OpenMS::SiriusMzTabWriter::SiriusAdapterIdentification":
        SiriusAdapterIdentification() nogil except +
        SiriusAdapterIdentification(SiriusAdapterIdentification) nogil except + # wrap-ignore

        double mz
        double rt
        String native_id
        int scan_index
        int scan_number
        String feature_id
        libcpp_vector[SiriusAdapterHit] hits

    cdef cppclass SiriusAdapterHit "OpenMS::SiriusMzTabWriter::SiriusAdapterHit":
        SiriusAdapterHit() nogil except +
        SiriusAdapterHit(SiriusAdapterHit) nogil except + # wrap-ignore

        String formula
        String adduct
        int rank
        double score
        double treescore
        double isoscore
        int explainedpeaks
        double explainedintensity

