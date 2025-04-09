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
        #  The representation of extra string data attached to a spectrum or chromatogram.

        StringDataArray() except + nogil 
        StringDataArray(StringDataArray &) except + nogil  # compiler

        Size size() except + nogil 
        void resize(size_t n) except + nogil 
        String& operator[](size_t) except + nogil  # wrap-ignore
        void clear() except + nogil 
        void push_back(String) except + nogil 
