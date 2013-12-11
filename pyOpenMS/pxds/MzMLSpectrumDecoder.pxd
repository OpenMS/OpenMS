from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *
from InterfaceDataStructures cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>" namespace "OpenMS":
    
    cdef cppclass MzMLSpectrumDecoder "OpenMS::MzMLSpectrumDecoder":
        MzMLSpectrumDecoder() nogil except +
        MzMLSpectrumDecoder(MzMLSpectrumDecoder) nogil except +
        void domParseChromatogram(libcpp_string in_, shared_ptr[Chromatogram] & cptr) nogil except +
        void domParseSpectrum(libcpp_string in_, shared_ptr[Spectrum] & cptr) nogil except +

