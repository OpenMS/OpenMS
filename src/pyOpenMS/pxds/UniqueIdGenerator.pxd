from String cimport *
from Types cimport *
from DateTime cimport *
from Param cimport *

cdef extern from "<OpenMS/CONCEPT/UniqueIdGenerator.h>" namespace "OpenMS":

    cdef cppclass UniqueIdGenerator:

        # protected
        UniqueIdGenerator() nogil except + # wrap-ignore
        # private
        UniqueIdGenerator(UniqueIdGenerator &) nogil except + # wrap-ignore

        # void setSeed(DateTime & time)

        UInt64 getUniqueId() nogil except +

        # Initializes random generator using the given DateTime instead of DateTime::now().  This is intended for debugging and testing.
        void setSeed(UInt64) nogil except +

        # Returns the seed
        UInt64 getSeed() nogil except +
