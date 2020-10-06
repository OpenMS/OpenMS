from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from MzTab cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMzTabWriter "OpenMS::SiriusMzTabWriter":
        SiriusMzTabWriter() nogil except + 
        SiriusMzTabWriter(SiriusMzTabWriter) nogil except + #wrap-ignore
        int extractScanIndex(const String& path) nogil except +
        int extractScanNumber(const String& path) nogil except +
        String extractFeatureId(const String& path) nogil except +
        void read(libcpp_vector[String] sirius_output_paths, String original_input_mzml, Size top_n_hits, MzTab & result) nogil except +
