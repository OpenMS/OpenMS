from Types cimport *
from String cimport *

from streampos cimport *


cdef extern from "<OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>" namespace "OpenMS":

    cdef cppclass IndexedMzMLDecoder:
        # wrap-doc:
            #   A class to analyze indexedmzML files and extract the offsets of individual tags
            #   -----
            #   Specifically, this class allows one to extract the offsets of the <indexList>
            #   tag and of all <spectrum> and <chromatogram> tag using the indices found at
            #   the end of the indexedmzML XML structure
            #   -----
            #   While findIndexListOffset tries extracts the offset of the indexList tag from
            #   the last 1024 bytes of the file, this offset allows the function parseOffsets
            #   to extract all elements contained in the <indexList> tag and thus get access
            #   to all spectra and chromatogram offsets

        IndexedMzMLDecoder() nogil except +
        IndexedMzMLDecoder(IndexedMzMLDecoder &) nogil except +
        int parseOffsets(String in_, int indexoffset, 
                libcpp_vector[ libcpp_pair[ libcpp_string, streampos] ]& spectra_offsets,
                libcpp_vector[ libcpp_pair[ libcpp_string, streampos] ]& chromatograms_offsets) nogil except + #wrap-ignore
        streampos findIndexListOffset(String in_, int buffersize) nogil except +
            # wrap-doc:
                #   Tries to extract the indexList offset from an indexedmzML
                #   -----
                #   This function reads by default the last few (1024) bytes of the given
                #   input file and tries to read the content of the <indexListOffset> tag
                #   The idea is that somewhere in the last parts of the file specified by the
                #   input string, the string <indexListOffset>xxx</indexListOffset> occurs
                #   This function returns the xxx part converted to an integer
                #   -----
                #   Since this function cannot determine where it will start reading
                #   the XML, no regular XML parser can be used for this. Therefore it uses
                #   regex to do its job. It matches the <indexListOffset> part and any
                #   numerical characters that follow
                #   -----
                #   :param in: Filename of the input indexedmzML file
                #   :param buffersize: How many bytes of the input file should be searched for the tag
                #   :returns: A positive integer containing the content of the indexListOffset tag, returns -1 in case of failure no tag was found (you can re-try with a larger buffersize but most likely its not an indexed mzML). Using -1 is what the reference docu recommends: http://en.cppreference.com/w/cpp/io/streamoff
                #   :raises:
                #     Exception: FileNotFound is thrown if file cannot be found
                #   :raises:
                #     Exception: ParseError if offset cannot be parsed
