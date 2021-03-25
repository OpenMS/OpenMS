from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from StringList cimport *
from MzTab cimport *
from CsvFile cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMzTabWriter "OpenMS::SiriusMzTabWriter":
        SiriusMzTabWriter() nogil except + 
        SiriusMzTabWriter(SiriusMzTabWriter) nogil except + #wrap-ignore

        int extractScanIndex(const String& path) nogil except +
        int extractScanNumber(const String& path) nogil except +
        String extractFeatureId(const String& path) nogil except +
        libcpp_map[ String, Size ] extract_columnname_to_columnindex(CsvFile& csvfile) nogil except +

        SiriusMzTabWriter_SiriusSpectrumMSInfo extractSpectrumMSInfo(const String& single_sirius_path) nogil except +

        void read(libcpp_vector[String] sirius_output_paths, String original_input_mzml, Size top_n_hits, MzTab & result) nogil except +

    cdef cppclass SiriusMzTabWriter_SiriusSpectrumMSInfo "OpenMS::SiriusMzTabWriter::SiriusSpectrumMSInfo":

      SiriusMzTabWriter_SiriusSpectrumMSInfo() nogil except +
      SiriusMzTabWriter_SiriusSpectrumMSInfo(SiriusMzTabWriter_SiriusSpectrumMSInfo) nogil except +

      StringList ext_n_id
      double ext_mz
      double ext_rt