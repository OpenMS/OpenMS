from String cimport *
from Types cimport *
from DateTime cimport *
from Param cimport *

cdef extern from "<OpenMS/CONCEPT/UniqueIdGenerator.h>" namespace "OpenMS":

    cdef cppclass UniqueIdGenerator:

        # protected
        UniqueIdGenerator() except + nogil  # wrap-ignore
        # private
        UniqueIdGenerator(UniqueIdGenerator &) except + nogil  # wrap-ignore

        # void setSeed(DateTime & time)

        UInt64 getUniqueId() except + nogil 

        # Initializes random generator using the given DateTime instead of DateTime::now().  This is intended for debugging and testing.
        void setSeed(UInt64) except + nogil 

        # Returns the seed
        UInt64 getSeed() except + nogil 
