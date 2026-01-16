from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoDescription cimport *

# see ../addons/IntegerDataArray.pyx

cdef extern from "<OpenMS/METADATA/DataArrays.h>" namespace "OpenMS::DataArrays":

    cdef cppclass IntegerDataArray(MetaInfoDescription):
        # wrap-inherits:
        #  MetaInfoDescription
        #
        # wrap-doc:
        #  The representation of extra integer data attached to a spectrum or chromatogram.
        #  Raw data access is provided by `get_peaks` and `set_peaks`, which yields numpy arrays.
        #  Used for storing per-peak integer annotations.

        IntegerDataArray() except + nogil  # wrap-doc:Default constructor
        IntegerDataArray(IntegerDataArray &) except + nogil  # wrap-doc:Copy constructor

        bool operator==(IntegerDataArray) except + nogil  # wrap-doc:Equality operator
        bool operator!=(IntegerDataArray) except + nogil  # wrap-doc:Inequality operator

        Size size() except + nogil  # wrap-doc:Returns the number of elements in the array
        void resize(size_t n) except + nogil  # wrap-doc:Resizes the array to contain n elements
        void reserve(size_t n) except + nogil  # wrap-doc:Reserves space for n elements to avoid repeated memory allocation
        Int& operator[](size_t) except + nogil  # wrap-ignore
        void clear() except + nogil  # wrap-doc:Removes all elements from the array
        void push_back(Int) except + nogil  # wrap-doc:Adds an integer value to the end of the array

        libcpp_vector[int].iterator begin() except + nogil   # wrap-ignore
        libcpp_vector[int].iterator end()   except + nogil   # wrap-ignore


