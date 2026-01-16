from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoDescription cimport *

# see ../addons/FloatDataArray.pyx

cdef extern from "<OpenMS/METADATA/DataArrays.h>" namespace "OpenMS::DataArrays":

    cdef cppclass FloatDataArray(MetaInfoDescription):
        # wrap-inherits:
        #  MetaInfoDescription

        # wrap-doc:
        #  The representation of extra float data attached to a spectrum or chromatogram.
        #  Raw data access is provided by `get_peaks` and `set_peaks`, which yields numpy arrays.
        #  Commonly used for storing ion mobility values or other per-peak float annotations.

        FloatDataArray() except + nogil  # wrap-doc:Default constructor
        FloatDataArray(FloatDataArray &) except + nogil  # wrap-doc:Copy constructor

        bool operator==(FloatDataArray) except + nogil  # wrap-doc:Equality operator
        bool operator!=(FloatDataArray) except + nogil  # wrap-doc:Inequality operator

        Size size() except + nogil  # wrap-doc:Returns the number of elements in the array
        void resize(size_t n) except + nogil  # wrap-doc:Resizes the array to contain n elements
        void reserve(size_t n) except + nogil  # wrap-doc:Reserves space for n elements to avoid repeated memory allocation
        float& operator[](size_t) except + nogil  # wrap-ignore
        void clear() except + nogil  # wrap-doc:Removes all elements from the array
        void push_back(float) except + nogil  # wrap-doc:Adds a float value to the end of the array

        libcpp_vector[float].iterator begin() nogil # wrap-ignore
        libcpp_vector[float].iterator end()   nogil # wrap-ignore
        void assign(float*, float*) # wrap-ignore


