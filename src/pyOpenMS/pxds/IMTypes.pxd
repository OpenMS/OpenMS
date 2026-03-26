from Types cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from libcpp.string cimport string as libcpp_utf8_string
from libcpp.string cimport string as libcpp_utf8_output_string

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
        CCS,
        SIZE_OF_DRIFTTIMEUNIT

    cdef enum IMFormat:
        NONE,
        IM_PEAK,
        IM_SPECTRUM,
        MIXED,
        UNKNOWN,
        SIZE_OF_IMFORMAT

    cdef enum IMPeakType:
        IM_PROFILE,
        IM_CENTROIDED,
        UNKNOWN,
        SIZE_OF_IMPEAKTYPE

# COMMENT: wrap static functions
cdef extern from "<OpenMS/IONMOBILITY/IMTypes.h>" namespace "OpenMS":

    # static fxn
    DriftTimeUnit toDriftTimeUnit(const libcpp_utf8_string& dtu_string) except + nogil  # wrap-attach:IMTypes
    libcpp_utf8_output_string driftTimeUnitToString(const DriftTimeUnit value) except + nogil  # wrap-attach:IMTypes

    IMFormat toIMFormat(const libcpp_utf8_string& IM_format) except + nogil  # wrap-attach:IMTypes
    libcpp_utf8_output_string imFormatToString(const IMFormat value) except + nogil  # wrap-attach:IMTypes

    IMPeakType toIMPeakType(const libcpp_utf8_string& im_peak_type) except + nogil  # wrap-attach:IMTypes
    libcpp_utf8_output_string imPeakTypeToString(const IMPeakType value) except + nogil  # wrap-attach:IMTypes

