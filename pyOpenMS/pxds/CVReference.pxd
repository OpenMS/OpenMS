from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVReference.h>" namespace "OpenMS":
    
    cdef cppclass CVReference "OpenMS::CVReference":
        CVReference() nogil except +
        CVReference(CVReference) nogil except +
        void setName(String &name) nogil except +
        String  getName() nogil except +
        void setIdentifier(String &identifier) nogil except +
        String  getIdentifier() nogil except +
        bool operator==(CVReference &rhs) nogil except +
        bool operator!=(CVReference &rhs) nogil except +

