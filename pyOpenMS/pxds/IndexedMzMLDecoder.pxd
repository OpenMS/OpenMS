from Types cimport *
from String cimport *

ctypedef libcpp_vector[ libcpp_pair[ libcpp_string, long] ] OffsetVector

cdef extern from "<OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>" namespace "OpenMS":
    
    cdef cppclass IndexedMzMLDecoder "OpenMS::IndexedMzMLDecoder":
        IndexedMzMLDecoder(IndexedMzMLDecoder) nogil except + #wrap-ignore
        int parseOffsets(String in_, int indexoffset, 
                         libcpp_vector[ libcpp_pair[ libcpp_string, long] ] & spectra_offsets,
                         libcpp_vector[ libcpp_pair[ libcpp_string, long] ] & chromatograms_offsets) nogil except +
        int findIndexListOffset(String in_, int buffersize) nogil except +

