from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/SpectrumNativeIDParser.h>" namespace "OpenMS":

    cdef cppclass SpectrumNativeIDParser:
        # wrap-doc:
        #  Parser for extracting scan numbers from spectrum native IDs

        # SpectrumNativeIDParser is a utility class with only static methods
        # No constructor needed for Python binding since we only need static methods

        pass

# Static methods exposed at module level
cdef extern from "<OpenMS/METADATA/SpectrumNativeIDParser.h>" namespace "OpenMS::SpectrumNativeIDParser":

    Int extractScanNumber(const String& native_id, const String& native_id_type_accession) except + nogil
    # wrap-attach:
    #   SpectrumNativeIDParser
    # wrap-doc:
    #   Extract the scan number from the native ID using a CV accession
    #
    #   :param native_id: Spectrum native ID string
    #   :param native_id_type_accession: CV accession specifying the native ID format (e.g., "MS:1000768" for Thermo, "MS:1000770" for WIFF)
    #   :returns: Scan number of the spectrum (or -1 on failure to extract)
    #
    #   Supported CV accessions:
    #     - MS:1000768, MS:1000769, MS:1000771, MS:1000772, MS:1000776: scan=NUMBER format
    #     - MS:1000770: WIFF format (returns cycle * 1000 + experiment)
    #     - MS:1000773, MS:1000775: file=NUMBER format
    #     - MS:1000774: index=NUMBER format (returns index + 1)
    #     - MS:1001508: scanId=NUMBER format
    #     - MS:1000777: spectrum=NUMBER format
    #     - MS:1001530: plain NUMBER format

    String getRegExFromNativeID "OpenMS::SpectrumNativeIDParser::getRegExFromNativeID" (const String& native_id) except + nogil
    # wrap-attach:
    #   SpectrumNativeIDParser
    # wrap-doc:
    #   Determine the regular expression to extract scan/index numbers from native IDs
    #
    #   :param native_id: A native ID string to analyze
    #   :returns: Regular expression string with named group that matches the scan or index number

    bool isNativeID "OpenMS::SpectrumNativeIDParser::isNativeID" (const String& id) except + nogil
    # wrap-attach:
    #   SpectrumNativeIDParser
    # wrap-doc:
    #   Check if a spectrum identifier is a native ID from a vendor file
    #
    #   :param id: Spectrum identifier string to check
    #   :returns: True if the string matches a known native ID prefix pattern
    #
    #   Recognized prefixes: scan=, scanId=, scanID=, controllerType=, function=, sample=, index=, spectrum=, file=
