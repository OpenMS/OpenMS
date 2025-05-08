from Types cimport *
from Matrix cimport *

cdef extern from "<OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>" namespace "OpenMS":
    
    cdef cppclass NonNegativeLeastSquaresSolver "OpenMS::NonNegativeLeastSquaresSolver":
        NonNegativeLeastSquaresSolver() except + nogil 
        NonNegativeLeastSquaresSolver(NonNegativeLeastSquaresSolver &) except + nogil 
        Int solve(Matrix[ double ] & A, Matrix[ double ] & b, Matrix[ double ] & x) except + nogil 

cdef extern from "<OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>" namespace "OpenMS::NonNegativeLeastSquaresSolver":
    cdef enum RETURN_STATUS "OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS":
        #wrap-attach:
        #   NonNegativeLeastSquaresSolver
        SOLVED
        ITERATION_EXCEEDED

