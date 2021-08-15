from Types cimport *
from Matrix cimport *

cdef extern from "<OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>" namespace "OpenMS":
    
    cdef cppclass NonNegativeLeastSquaresSolver "OpenMS::NonNegativeLeastSquaresSolver":
        NonNegativeLeastSquaresSolver() nogil except +
        NonNegativeLeastSquaresSolver(NonNegativeLeastSquaresSolver &) nogil except +
        Int solve(Matrix[ double ] & A, Matrix[ double ] & b, Matrix[ double ] & x) nogil except +

cdef extern from "<OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>" namespace "OpenMS::NonNegativeLeastSquaresSolver":
    cdef enum RETURN_STATUS "OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS":
        #wrap-attach:
        #    NonNegativeLeastSquaresSolver
        SOLVED
        ITERATION_EXCEEDED

