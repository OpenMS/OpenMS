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

        # wrap-doc:
        #  The representation of extra string data attached to a spectrum or chromatogram.
        #  Commonly used for storing ion annotation names or other per-peak string annotations.

        StringDataArray() except + nogil  # wrap-doc:Default constructor
        StringDataArray(StringDataArray &) except + nogil  # wrap-doc:Copy constructor

        bool operator==(StringDataArray) except + nogil  # wrap-doc:Equality operator
        bool operator!=(StringDataArray) except + nogil  # wrap-doc:Inequality operator

        Size size() except + nogil  # wrap-doc:Returns the number of elements in the array
        void resize(size_t n) except + nogil  # wrap-doc:Resizes the array to contain n elements
        # Implemented in StringDataArray.pyx
        String& operator[](size_t) except + nogil  # wrap-ignore
        void clear() except + nogil  # wrap-doc:Removes all elements from the array
        void push_back(String) except + nogil  # wrap-doc:Adds a string value to the end of the array
