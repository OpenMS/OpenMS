from DataValue cimport *

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS":

    cdef cppclass CVTerm:
         CVTerm()
         CVTerm(CVTerm)
         DataValue getValue()


