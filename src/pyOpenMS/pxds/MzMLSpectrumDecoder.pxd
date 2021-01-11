from Types cimport *
from String cimport *
from InterfaceDataStructures cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>" namespace "OpenMS":
    
    cdef cppclass MzMLSpectrumDecoder "OpenMS::MzMLSpectrumDecoder":
        MzMLSpectrumDecoder() nogil except +
        MzMLSpectrumDecoder(MzMLSpectrumDecoder) nogil except +
        void domParseChromatogram(String in_, shared_ptr[Chromatogram] & cptr) nogil except +
        void domParseSpectrum(String in_, shared_ptr[Spectrum] & cptr) nogil except +
        void setSkipXMLChecks(bool only) nogil except +

