from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVReference.h>" namespace "OpenMS":
    
    cdef cppclass CVReference "OpenMS::CVReference":
        CVReference() nogil except +
        CVReference(CVReference &) nogil except +
        void setName(const String &name) nogil except + # wrap-doc:Sets the name of the CV reference
        String getName() nogil except + # wrap-doc:Returns the name of the CV reference
        void setIdentifier(const String &identifier) nogil except + # wrap-doc:Sets the CV identifier which is referenced
        String getIdentifier() nogil except + # wrap-doc:Returns the CV identifier which is referenced
        bool operator==(CVReference &rhs) nogil except +
        bool operator!=(CVReference &rhs) nogil except +

