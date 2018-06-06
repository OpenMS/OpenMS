from String cimport *
from Types cimport *

cdef extern from "<OpenMS/CONCEPT/UniqueIdInterface.h>" namespace "OpenMS":

    cdef cppclass UniqueIdInterface:
        # wrap-ignore
        # no-pxd-import

        UniqueIdInterface() nogil except +
        UniqueIdInterface(UniqueIdInterface) nogil except +

        Size getUniqueId() nogil except +
        Size clearUniqueId() nogil except +
        void swap(UniqueIdInterface) nogil except + # wrap-ignore
        Size hasValidUniqueId() nogil except +
        Size hasInvalidUniqueId() nogil except +
        void setUniqueId(UInt64 rhs) nogil except +
        Size ensureUniqueId() nogil except +
        
        # overloading Cython Issue
        # TODO: Mismatch between C++ return type ([u'Size']) and Python return type (['void']) in function public setUniqueId:
        # Size setUniqueId() nogil except +
        # void setUniqueId(const String & rhs) nogil except +

        bool isValid(UInt64 unique_id) nogil except +

cdef extern from "<OpenMS/CONCEPT/UniqueIdInterface.h>" namespace "OpenMS::UniqueIdInterface":

        Size setUniqueId() # wrap-ignore




