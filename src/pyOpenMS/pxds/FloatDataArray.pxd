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
        #   The representation of extra float data attached to a spectrum or chromatogram.
        #   Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays

        FloatDataArray() nogil except +
        FloatDataArray(FloatDataArray) nogil except + #wrap-ignore

        Size size() nogil except +
        void resize(size_t n) nogil except +
        void reserve(size_t n) nogil except +
        float& operator[](int) nogil except + # wrap-ignore
        void clear() nogil except +
        void push_back(float) nogil except +

        libcpp_vector[float].iterator begin() nogil # wrap-ignore
        libcpp_vector[float].iterator end()   nogil # wrap-ignore
        void assign(float*, float*) # wrap-ignore


