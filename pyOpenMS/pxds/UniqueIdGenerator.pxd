from String cimport *
from Types cimport *
from DateTime cimport *
from Param cimport *

cdef extern from "<OpenMS/CONCEPT/UniqueIdGenerator.h>" namespace "OpenMS":

    cdef cppclass UniqueIdGenerator:

        # void setSeed(DateTime & time)

        UInt64 getUniqueId() nogil except +

        # Initializes random generator using the given DateTime instead of DateTime::now().  This is intended for debugging and testing.
        void setSeed(DateTime) nogil except +

        # Returns a summary of internal status
        Param getInfo() nogil except +

