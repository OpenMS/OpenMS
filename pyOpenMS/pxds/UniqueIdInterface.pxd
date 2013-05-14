from String cimport *
from Types cimport *

cdef extern from "<OpenMS/CONCEPT/UniqueIdInterface.h>" namespace "OpenMS":

    cdef cppclass UniqueIdInterface:
        # wrap-ignore

        UniqueIdInterface()
        UniqueIdInterface(UniqueIdInterface)

        Size getUniqueId() nogil except +
        Size clearUniqueId() nogil except +
        void swap(UniqueIdInterface) nogil except + # wrap-ignore
        Size hasValidUniqueId() nogil except +
        Size hasInvalidUniqueId() nogil except +
        void setUniqueId(UInt64 rhs) nogil except +
        Size ensureUniqueId() nogil except +
        
        # overloading Cython Issue
        # Size setUniqueId() nogil except +
        # void setUniqueId(String & rhs) nogil except +

        bool isValid(UInt64 unique_id)

cdef extern from "<OpenMS/CONCEPT/UniqueIdInterface.h>" namespace "OpenMS::UniqueIdInterface":

        Size setUniqueId()  # wrap-ignore




