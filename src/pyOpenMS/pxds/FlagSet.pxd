from libcpp cimport bool
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/FlagSet.h>" namespace "OpenMS":
    
    cdef cppclass FlagSet[ENUM]:
        # wrap-ignore
        # wrap-doc:
        #  Stores and handles combinations of enum values
        #  
        #  This class stores combinations of enum values as bits in an integer.
        #  This is used for example in QCBase::Status to store requirements.
        
        FlagSet() except + nogil 
        FlagSet(ENUM en) except + nogil 
        FlagSet(FlagSet[ENUM] &) except + nogil 
        
        bool operator==(FlagSet[ENUM] & stat) except + nogil 
        
        FlagSet[ENUM] operator&(ENUM en) except + nogil 
        FlagSet[ENUM] operator&(FlagSet[ENUM] & rhs) except + nogil 
        
        FlagSet[ENUM] operator|(ENUM en) except + nogil 
        FlagSet[ENUM] operator|(FlagSet[ENUM] & rhs) except + nogil 
        
        bool isSuperSetOf(FlagSet[ENUM] & required) except + nogil 
