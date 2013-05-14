from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DistanceMatrix.h>" namespace "OpenMS":
    
    cdef cppclass DistanceMatrix[Value]:
        # wrap-instances:
        #   DistanceMatrix := DistanceMatrix[float]
        DistanceMatrix() nogil except +
        DistanceMatrix(DistanceMatrix) nogil except +
        DistanceMatrix(size_t dimensionsize, Value value) nogil except +
        # ValueType operator()(size_t i, size_t j) nogil except +
        # ValueType operator()(size_t i, size_t j) nogil except +
        Value getValue(size_t i, size_t j) nogil except +
        void setValue(size_t i, size_t j, Value value) nogil except +
        void setValueQuick(size_t i, size_t j, Value value) nogil except +
        void clear() nogil except +
        void resize(size_t dimensionsize, Value value) nogil except +
        void reduce(size_t j) nogil except +
        size_t dimensionsize() nogil except +
        void updateMinElement() nogil except +
        bool operator==(DistanceMatrix[ Value ] &rhs) nogil except +
        libcpp_pair[ size_t, size_t ] getMinElementCoordinates() nogil except +
        # TEMPLATE # std::ostream  operator[[(std::ostream &os, DistanceMatrix[ Value ] &matrix) nogil except +

