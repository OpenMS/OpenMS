from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/LinearInterpolation.h>" namespace "OpenMS::Math":
    
    cdef cppclass LinearInterpolation[KeyType,ValueType]:
        # wrap-instances:
        #   LinearInterpolation := LinearInterpolation[double, double]
        LinearInterpolation() nogil except +
        LinearInterpolation(LinearInterpolation) nogil except +
        ValueType value(KeyType arg_pos) nogil except +
        void addValue(KeyType arg_pos, ValueType arg_value) nogil except +
        ValueType derivative(KeyType arg_pos) nogil except +
        # TODO does this work ?
        # libcpp_vector[ValueType]  getData() nogil except +
        # void setData(libcpp_vector[ValueType] & data) nogil except +
        bool empty() nogil except +
        KeyType key2index(KeyType pos) nogil except +
        KeyType index2key(KeyType pos) nogil except +
        KeyType  getScale() nogil except +
        void setScale(KeyType & scale) nogil except +
        KeyType  getOffset() nogil except +
        void setOffset(KeyType & offset) nogil except +
        void setMapping(KeyType & scale, KeyType & inside, KeyType & outside) nogil except +
        void setMapping(KeyType & inside_low, KeyType & outside_low, KeyType & inside_high, KeyType & outside_high) nogil except +
        KeyType  getInsideReferencePoint() nogil except +
        KeyType  getOutsideReferencePoint() nogil except +
        KeyType supportMin() nogil except +
        KeyType supportMax() nogil except +
        LinearInterpolation(KeyType scale, KeyType offset) nogil except +

