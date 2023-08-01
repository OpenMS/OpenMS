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
        # const_reference operator()(size_t i, size_t j) except + nogil 
        # reference operator()(size_t i, size_t j) except + nogil 
        # const_reference getValue(size_t i, size_t j) except + nogil 
        ValueT getValue(size_t i, size_t j) nogil
        void setValue(size_t i, size_t j, ValueT value) nogil
        ## The following two lines introduce an odd bug:
        # static PyObject *__pyx_convert_vector_to_py_double is declared twice by Cython:
        # TODO look into Cython Bug
        # libcpp_vector[ValueT] row(size_t i) except + nogil 
        # libcpp_vector[ValueT] col(size_t i) except + nogil 
        libcpp_vector[ValueT] asVector() except + nogil  # wrap-ignore
        void clear() except + nogil 
        void resize(size_t i, size_t j, ValueT value) except + nogil 
        void resize(libcpp_pair[ size_t, size_t ] & size_pair, ValueT value) except + nogil 
        size_t rows() nogil
        size_t cols() nogil
        libcpp_pair[ size_t, size_t ] sizePair() except + nogil 
        size_t index(size_t row, size_t col) except + nogil 
        libcpp_pair[ size_t, size_t ] indexPair(size_t index) except + nogil 
        size_t colIndex(size_t index) except + nogil  # wrap-doc:Calculate the column from an index into the underlying vector. Note that Matrix uses the (row,column) lexicographic ordering for indexing
        size_t rowIndex(size_t index) except + nogil  # wrap-doc:Calculate the row from an index into the underlying vector. Note that Matrix uses the (row,column) lexicographic ordering for indexing
        ## bool operator==(Matrix & rhs) except + nogil 
        ## bool operator<(Matrix & rhs) except + nogil 
        # TEMPLATE # void setMatrix(ValueType matrix) except + nogil 
        # TEMPLATE # # NAMESPACE # std::ostream  operator<[(std::ostream & os, Matrix[ Value ] & matrix) except + nogil 
        #  MatrixUnsignedInt := Matrix[unsigned int]

