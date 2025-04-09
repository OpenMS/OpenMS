from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/DATASTRUCTURES/Matrix.h>" namespace "OpenMS":
    
    cdef cppclass Matrix[ValueT]:
        # wrap-instances:
        #  MatrixDouble := Matrix[double]

        Matrix() except + nogil 
        Matrix(Matrix[ValueT]) except + nogil 
        Matrix(size_t rows, size_t cols, ValueT value) except + nogil 
        ValueT getValue(size_t i, size_t j) nogil
        void setValue(size_t i, size_t j, ValueT value) nogil
        size_t rows() nogil
        size_t cols() nogil
        size_t size() nogil
        int innerStride() nogil # wrap-ignore
        int outerStride() nogil # wrap-ignore
        bool rowMajor() nogil # wrap-ignore
        void resize(size_t rows, size_t cols) nogil
        ValueT* data() nogil # wrap-ignore
        ## bool operator==(Matrix & rhs) except + nogil 
        ## bool operator<(Matrix & rhs) except + nogil 
        # TEMPLATE # void setMatrix(ValueType matrix) except + nogil 
        # TEMPLATE # # NAMESPACE # std::ostream  operator<[(std::ostream & os, Matrix[ Value ] & matrix) except + nogil 
        #  MatrixUnsignedInt := Matrix[unsigned int]

