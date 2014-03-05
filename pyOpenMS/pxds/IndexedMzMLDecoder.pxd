from Types cimport *
from String cimport *

from streampos cimport *


cdef extern from "<OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>" namespace "OpenMS":

    cdef cppclass IndexedMzMLDecoder:
        IndexedMzMLDecoder() nogil except +
        IndexedMzMLDecoder(IndexedMzMLDecoder) nogil except +
        int parseOffsets(String in_, int indexoffset, 
                libcpp_vector[ libcpp_pair[ libcpp_string, streampos] ]& spectra_offsets,
                libcpp_vector[ libcpp_pair[ libcpp_string, streampos] ]& chromatograms_offsets) nogil except + #wrap-ignore
        streampos findIndexListOffset(String in_, int buffersize) nogil except +
