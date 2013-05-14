from String cimport *
from Types cimport *
from DateTime cimport *

cdef extern from "<OpenMS/CONCEPT/UniqueIdGenerator.h>" namespace "OpenMS":

    cdef cppclass UniqueIdGenerator:

        void setSeed(DateTime & time)




