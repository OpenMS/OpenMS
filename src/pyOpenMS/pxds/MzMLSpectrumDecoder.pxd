from Types cimport *
from String cimport *
from InterfaceDataStructures cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>" namespace "OpenMS":
    
    cdef cppclass MzMLSpectrumDecoder "OpenMS::MzMLSpectrumDecoder":
        # wrap-doc:
                #   A class to decode input strings that contain an mzML chromatogram or spectrum tag
                #   -----
                #   It uses xercesc to parse a string containing either a exactly one mzML
                #   spectrum or chromatogram (from <chromatogram> to </chromatogram> or
                #   <spectrum> to </spectrum> tag). It returns the data contained in the
                #   binaryDataArray for Intensity / mass-to-charge or Intensity / time
                
        MzMLSpectrumDecoder() nogil except + # compiler
        MzMLSpectrumDecoder(MzMLSpectrumDecoder &) nogil except + # compiler
        void domParseChromatogram(String in_, shared_ptr[Chromatogram] & cptr) nogil except +
            # wrap-doc:
                #   Extract data from a string which contains a full mzML chromatogram
                #   -----
                #   Extracts data from the input string which is expected to contain exactly
                #   one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
                #   function will extract the contained binaryDataArray and provide the
                #   result as Chromatogram
                #   -----
                #   :param in: Input string containing the raw XML
                #   :param cptr: Resulting chromatogram

        void domParseSpectrum(String in_, shared_ptr[Spectrum] & cptr) nogil except +
            # wrap-doc:
                #   Extract data from a string which contains a full mzML spectrum
                #   -----
                #   Extracts data from the input string which is expected to contain exactly
                #   one <spectrum> tag (from <spectrum> to </spectrum>). This function will
                #   extract the contained binaryDataArray and provide the result as Spectrum
                #   -----
                #   :param in: Input string containing the raw XML
                #   :param cptr: Resulting spectrum

        void setSkipXMLChecks(bool only) nogil except + # wrap-doc:Whether to skip some XML checks (e.g. removing whitespace inside base64 arrays) and be fast instead

