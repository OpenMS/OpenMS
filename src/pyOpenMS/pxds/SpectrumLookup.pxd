from Types cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/METADATA/SpectrumLookup.h>" namespace "OpenMS":

    cdef cppclass SpectrumLookup:

        SpectrumLookup() nogil except +
        # SpectrumLookup(SpectrumLookup) nogil except + # private

        # TODO: Found attribute in C++ but not in pxd: variable public default_scan_regexp
        # TODO: Found attribute in C++ but not in pxd: variable public reference_formats
        # TODO: Found attribute in C++ but not in pxd: variable public rt_tolerance
        double rt_tolerance

        bool empty() nogil except +

        void readSpectra(MSExperiment spectra, String scan_regexp) nogil except +

        Size findByRT(double rt) nogil except +

        Size findByNativeID(String native_id) nogil except +
        
        # This is a base class: we cannot overload methods since Cython doesn't
        # allow inheritance of overloaded methods.

        # Size findByIndex(Size index) nogil except +

        Size findByIndex(Size index, bool count_from_one) nogil except +

        Size findByScanNumber(Size scan_number) nogil except +

        Size findByReference(String spectrum_ref) nogil except +

        void addReferenceFormat(String regexp) nogil except +

        # # NAMESPACE # Int extractScanNumber(const String & native_id, boost::regex & scan_regexp, bool no_error) nogil except +

