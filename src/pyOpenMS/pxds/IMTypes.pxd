from Types cimport *
from MSExperiment cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/IONMOBILITY/IMTypes.h>" namespace "OpenMS":

    cdef cppclass IMTypes:

      IMFormat determineIMFormat(const MSExperiment& exp) nogil except +
      IMFormat determineIMFormat(const MSSpectrum& spec) nogil except +

cdef extern from "<OpenMS/IONMOBILITY/IMTypes.h>" namespace "OpenMS":

    cdef enum DriftTimeUnit:
        NONE,                      
        MILLISECOND,               
        VSSC,                      
        FAIMS_COMPENSATION_VOLTAGE,
        SIZE_OF_DRIFTTIMEUNIT

    cdef enum IMFormat:
        NONE,            
        CONCATENATED,    
        MULTIPLE_SPECTRA,
        MIXED,           
        SIZE_OF_IMFORMAT

# COMMENT: wrap static functions
cdef extern from "<OpenMS/IONMOBILITY/IMTypes.h>" namespace "OpenMS":
        
    # static fxn
    DriftTimeUnit toDriftTimeUnit(const libcpp_string& dtu_string) nogil except +

    # OPENMS_DLLAPI DriftTimeUnit toDriftTimeUnit(const std::string& dtu_string);
    # OPENMS_DLLAPI const std::string& toString(const DriftTimeUnit value);
    # OPENMS_DLLAPI extern const std::string NamesOfIMFormat[(size_t) IMFormat::SIZE_OF_IMFORMAT];
    # /// convert an entry in NamesOfIMFormat[] to IMFormat enum
    # /// @throws Exception::InvalidValue if @p IM_format is not contained in NamesOfIMFormat[]
    # OPENMS_DLLAPI IMFormat toIMFormat(const std::string& IM_format);
    # /// convert an IMFormat enum to String
    # /// @throws Exception::InvalidValue if @p value is SIZE_OF_IMFORMAT
    # OPENMS_DLLAPI const std::string& toString(const IMFormat value);


