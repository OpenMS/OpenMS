from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Map cimport *
from MzTab cimport *
from CsvFile cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMzTabWriter "OpenMS::SiriusMzTabWriter":
        SiriusMzTabWriter() nogil except + 
        SiriusMzTabWriter(SiriusMzTabWriter) nogil except +

        int extractScanIndex(const String& path) nogil except +
        int extractScanNumber(const String& path) nogil except +
        String extractFeatureId(const String& path) nogil except +
        # Map[ String, Size ] extract_columnname_to_columnindex(CsvFile& csvfile) nogil except +

        void read(libcpp_vector[String] sirius_output_paths, String original_input_mzml, Size top_n_hits, MzTab & result) nogil except +

cdef cppclass SiriusMzTabWriter_SiriusSpectrumMSInfo "OpenMS::SiriusMzTabWriter::SiriusSpectrumMSInfo":

   SiriusMzTabWriter_SiriusSpectrumMSInfo() nogil except +
   SiriusMzTabWriter_SiriusSpectrumMSInfo(SiriusMzTabWriter_SiriusSpectrumMSInfo) nogil except +