from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DistanceMatrix.h>" namespace "OpenMS":
    
    cdef cppclass DistanceMatrix[Value]:
        # wrap-instances:
        #  DistanceMatrix := DistanceMatrix[float]
        DistanceMatrix() except + nogil 
        DistanceMatrix(DistanceMatrix &) except + nogil 
        DistanceMatrix(size_t dimensionsize, Value value) except + nogil 

        # ValueType operator()(size_t i, size_t j) except + nogil 
        # ValueType operator()(size_t i, size_t j) except + nogil 
        Value getValue(size_t i, size_t j) except + nogil 
        void setValue(size_t i, size_t j, Value value) except + nogil 
        void setValueQuick(size_t i, size_t j, Value value) except + nogil 
        void clear() except + nogil 
        void resize(size_t dimensionsize, Value value) except + nogil 
        void reduce(size_t j) except + nogil 
        size_t dimensionsize() except + nogil 
        void updateMinElement() except + nogil 
        bool operator==(DistanceMatrix[ Value ] &rhs) except + nogil 
        libcpp_pair[ size_t, size_t ] getMinElementCoordinates() except + nogil 
        # TEMPLATE # std::ostream  operator[[(std::ostream &os, DistanceMatrix[ Value ] &matrix) except + nogil 

