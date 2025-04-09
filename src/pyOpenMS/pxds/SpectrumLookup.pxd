from Types cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/METADATA/SpectrumLookup.h>" namespace "OpenMS":

    cdef cppclass SpectrumLookup:

        SpectrumLookup() except + nogil 
        # private
        SpectrumLookup(SpectrumLookup &) except + nogil  # wrap-ignore

        double rt_tolerance

        bool empty() except + nogil  # wrap-doc:Check if any spectra were set

        void readSpectra(MSExperiment spectra, String scan_regexp) except + nogil 
        # wrap-doc:
                #  Read and index spectra for later look-up
                #  
                #  :param spectra: Container of spectra
                #  :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>". For example, "scan=(?<SCAN>\\d+)").

        Size findByRT(double rt) except + nogil 
        # wrap-doc:
                #  Look up spectrum by retention time (RT)
                #  
                #  :param rt: Retention time to look up
                #  :returns: Index of the spectrum that matched

        Size findByNativeID(String native_id) except + nogil 
        # wrap-doc:
                #  Look up spectrum by native ID
                #  
                #  :param native_id: Native ID to look up
                #  :returns: Index of the spectrum that matched
        
        # This is a base class: we cannot overload methods since Cython doesn't
        # allow inheritance of overloaded methods.

        # Size findByIndex(Size index) except + nogil 

        Size findByIndex(Size index, bool count_from_one) except + nogil 
         # wrap-doc:
                #  Look up spectrum by index (position in the vector of spectra)
                #  
                #  :param index: Index to look up
                #  :param count_from_one: Do indexes start counting at one (default zero)?
                #  :returns: Index of the spectrum that matched

        Size findByScanNumber(Size scan_number) except + nogil 
        # wrap-doc:
                #  Look up spectrum by scan number (extracted from the native ID)
                #  
                #  :param scan_number: Scan number to look up
                #  :returns: Index of the spectrum that matched

        Size findByReference(String spectrum_ref) except + nogil 
        # wrap-doc:
                #  Look up spectrum by reference
                #  
                #  :param spectrum_ref: Spectrum reference to parse
                #  :returns: Index of the spectrum that matched

        void addReferenceFormat(String regexp) except + nogil 
        # wrap-doc:
                #  Register a possible format for a spectrum reference
                #  
                #  :param regexp: Regular expression defining the format

        # # NAMESPACE # Int extractScanNumber(const String & native_id, boost::regex & scan_regexp, bool no_error) except + nogil 

        Int extractScanNumber(const String& native_id, const String& native_id_type_accession) except + nogil 
