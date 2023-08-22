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
        #
        # wrap-doc:
        #  The representation of extra float data attached to a spectrum or chromatogram.
        #  Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays

        FloatDataArray() except + nogil 
        FloatDataArray(FloatDataArray &) except + nogil  # compiler

        Size size() except + nogil 
        void resize(size_t n) except + nogil 
        void reserve(size_t n) except + nogil 
        float& operator[](size_t) except + nogil  # wrap-ignore
        void clear() except + nogil 
        void push_back(float) except + nogil 

        libcpp_vector[float].iterator begin() nogil # wrap-ignore
        libcpp_vector[float].iterator end()   nogil # wrap-ignore
        void assign(float*, float*) # wrap-ignore


