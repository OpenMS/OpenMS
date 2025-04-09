from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVReference.h>" namespace "OpenMS":
    
    cdef cppclass CVReference "OpenMS::CVReference":
        CVReference() except + nogil 
        CVReference(CVReference &) except + nogil 
        void setName(const String &name) except + nogil  # wrap-doc:Sets the name of the CV reference
        String getName() except + nogil  # wrap-doc:Returns the name of the CV reference
        void setIdentifier(const String &identifier) except + nogil  # wrap-doc:Sets the CV identifier which is referenced
        String getIdentifier() except + nogil  # wrap-doc:Returns the CV identifier which is referenced
        bool operator==(CVReference &rhs) except + nogil 
        bool operator!=(CVReference &rhs) except + nogil 

