from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoDescription cimport *

# see ../addons/StringDataArray.pyx

cdef extern from "<OpenMS/METADATA/DataArrays.h>" namespace "OpenMS::DataArrays":

    cdef cppclass StringDataArray(MetaInfoDescription):
        # wrap-inherits:
        #  MetaInfoDescription
        #
        # wrap-doc:
        #   The representation of extra string data attached to a spectrum or chromatogram.

        StringDataArray() nogil except +
        StringDataArray(StringDataArray) nogil except + #wrap-ignore

        Size size() nogil except +
        void resize(size_t n) nogil except +
        String& operator[](int) nogil except + # wrap-ignore
        void clear() nogil except +
        void push_back(String) nogil except +


