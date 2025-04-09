from String cimport *
from Types cimport *

cdef extern from "<OpenMS/CONCEPT/UniqueIdInterface.h>" namespace "OpenMS":

    cdef cppclass UniqueIdInterface:
        # wrap-ignore
        # no-pxd-import

        UniqueIdInterface() except + nogil  # wrap-doc:A base class defining a common interface for all classes having a unique id. Have a look at RichPeak2D for an example how to extend a class to support unique ids
        UniqueIdInterface(UniqueIdInterface &) except + nogil 

        Size getUniqueId() except + nogil  # wrap-doc:Returns the unique id
        Size clearUniqueId() except + nogil  # wrap-doc:Clear the unique id. The new unique id will be invalid.  Returns 1 if the unique id was changed, 0 otherwise
        void swap(UniqueIdInterface) except + nogil  # wrap-ignore
        Size hasValidUniqueId() except + nogil  # wrap-doc:Returns whether the unique id is valid.  Returns 1 if the unique id is valid, 0 otherwise
        Size hasInvalidUniqueId() except + nogil  # wrap-doc:Returns whether the unique id is invalid.  Returns 1 if the unique id is invalid, 0 otherwise
        void setUniqueId(UInt64 rhs) except + nogil  # wrap-doc:Assigns a new, valid unique id.  Always returns 1
        Size ensureUniqueId() except + nogil  # wrap-doc:Assigns a valid unique id, but only if the present one is invalid.  Returns 1 if the unique id was changed, 0 otherwise
        
        # overloading Cython Issue
        # TODO: Mismatch between C++ return type ([u'Size']) and Python return type (['void']) in function public setUniqueId:
        # Size setUniqueId() except + nogil 
        # void setUniqueId(const String & rhs) except + nogil 

        bool isValid(UInt64 unique_id) except + nogil  # wrap-doc:Returns true if the unique_id is valid, false otherwise

cdef extern from "<OpenMS/CONCEPT/UniqueIdInterface.h>" namespace "OpenMS::UniqueIdInterface":

        Size setUniqueId() # wrap-ignore
