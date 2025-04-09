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
        #  Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays

        IntegerDataArray() except + nogil 
        IntegerDataArray(IntegerDataArray &) except + nogil 

        Size size() except + nogil 
        void resize(size_t n) except + nogil 
        void reserve(size_t n) except + nogil 
        Int& operator[](size_t) except + nogil  # wrap-ignore
        void clear() except + nogil 
        void push_back(Int) except + nogil 

        libcpp_vector[int].iterator begin() except + nogil   # wrap-ignore
        libcpp_vector[int].iterator end()   except + nogil   # wrap-ignore


