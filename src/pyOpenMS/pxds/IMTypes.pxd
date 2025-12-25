from Types cimport *
from MSExperiment cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/IONMOBILITY/IMTypes.h>" namespace "OpenMS":

    cdef cppclass IMTypes:

      IMTypes() except + nogil 
      IMTypes(IMTypes &) except + nogil  # compiler
      IMFormat determineIMFormat(const MSExperiment& exp) except + nogil 
      IMFormat determineIMFormat(const MSSpectrum& spec) except + nogil 

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
    DriftTimeUnit toDriftTimeUnit(const libcpp_string& dtu_string) except + nogil  # wrap-attach:IMTypes
    libcpp_string toString(const DriftTimeUnit value) except + nogil  # wrap-attach:IMTypes

    IMFormat toIMFormat(const libcpp_string& IM_format) except + nogil  # wrap-attach:IMTypes
    libcpp_string toString(const IMFormat value) except + nogil  # wrap-attach:IMTypes


