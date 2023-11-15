from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from libcpp.map cimport map as libcpp_map
from String cimport *
from StringList cimport *
from MzTab cimport *
from CsvFile cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMzTabWriter:
        SiriusMzTabWriter() except + nogil  
        SiriusMzTabWriter(SiriusMzTabWriter &) except + nogil  # compiler

    cdef cppclass SiriusMzTabWriter_SiriusSpectrumMSInfo "OpenMS::SiriusMzTabWriter::SiriusSpectrumMSInfo":

      SiriusMzTabWriter_SiriusSpectrumMSInfo() except + nogil 
      SiriusMzTabWriter_SiriusSpectrumMSInfo(SiriusMzTabWriter_SiriusSpectrumMSInfo &) except + nogil  # compiler

      StringList ext_n_id
      double ext_mz
      double ext_rt

# wrap static method:
cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":

        int extractScanIndex(const String& path) except + nogil   # wrap-attach:SiriusMzTabWriter
        int extractScanNumber(const String& path) except + nogil  # wrap-attach:SiriusMzTabWriter
        String extractFeatureId(const String& path) except + nogil  # wrap-attach:SiriusMzTabWriter
        libcpp_map[ libcpp_string, Size ] extract_columnname_to_columnindex(CsvFile& csvfile) except + nogil  # wrap-attach:SiriusMzTabWriter

        SiriusMzTabWriter_SiriusSpectrumMSInfo extractSpectrumMSInfo(const String& single_sirius_path) except + nogil  # wrap-attach:SiriusMzTabWriter

        void read(libcpp_vector[ String ]& sirius_output_paths, 
                  const String& original_input_mzml, 
                  Size top_n_hits, 
                  MzTab& result) except + nogil  # wrap-attach:SiriusMzTabWriter
