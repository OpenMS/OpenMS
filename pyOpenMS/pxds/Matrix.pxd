from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/DATASTRUCTURES/Matrix.h>" namespace "OpenMS":
    
    cdef cppclass Matrix[ValueT]:
        # wrap-instances:
        #   MatrixDouble := Matrix[double]
        Matrix() nogil except +
        Matrix(Matrix[ValueT]) nogil except +
        Matrix(size_t rows, size_t cols, ValueT value) nogil except +
        # const_reference operator()(size_t i, size_t j) nogil except +
        # reference operator()(size_t i, size_t j) nogil except +
        # const_reference getValue(size_t i, size_t j) nogil except +
        ValueT getValue(size_t i, size_t j) nogil except +
        void setValue(size_t i, size_t j, ValueT value) nogil except +
        ## The following two lines introduce an odd bug:
        # static PyObject *__pyx_convert_vector_to_py_double is declared twice by Cython:
        # TODO look into Cython Bug
        # libcpp_vector[ValueT] row(size_t i) nogil except +
        # libcpp_vector[ValueT] col(size_t i) nogil except +
        void clear() nogil except +
        void resize(size_t i, size_t j, ValueT value) nogil except +
        void resize(libcpp_pair[ size_t, size_t ] & size_pair, ValueT value) nogil except +
        size_t rows() nogil except +
        size_t cols() nogil except +
        libcpp_pair[ size_t, size_t ] sizePair() nogil except +
        size_t index(size_t row, size_t col) nogil except +
        libcpp_pair[ size_t, size_t ] indexPair(size_t index) nogil except +
        size_t colIndex(size_t index) nogil except +
        size_t rowIndex(size_t index) nogil except +
        ## bool operator==(Matrix & rhs) nogil except +
        ## bool operator<(Matrix & rhs) nogil except +
        # TEMPLATE # void setMatrix(ValueType matrix) nogil except +
        # POINTER # gsl_matrix * toGslMatrix() nogil except +
        # TEMPLATE # # NAMESPACE # std::ostream  operator<[(std::ostream & os, Matrix[ Value ] & matrix) nogil except +
        #   MatrixUnsignedInt := Matrix[unsigned int]

