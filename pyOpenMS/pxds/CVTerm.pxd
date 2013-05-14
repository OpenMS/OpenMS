from DataValue cimport *

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS":

    cdef cppclass CVTerm:
         CVTerm()   nogil except +
         CVTerm(CVTerm)   nogil except +
         DataValue getValue()   nogil except +


