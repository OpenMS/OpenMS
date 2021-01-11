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
        #   The representation of extra integer data attached to a spectrum or chromatogram.
        #   Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays

        IntegerDataArray() nogil except +
        IntegerDataArray(IntegerDataArray) nogil except + #wrap-ignore

        Size size() nogil except +
        void resize(size_t n) nogil except +
        void reserve(size_t n) nogil except +
        Int& operator[](int) nogil except + # wrap-ignore
        void clear() nogil except +
        void push_back(Int) nogil except +

        libcpp_vector[int].iterator begin() nogil except +  # wrap-ignore
        libcpp_vector[int].iterator end()   nogil except +  # wrap-ignore


