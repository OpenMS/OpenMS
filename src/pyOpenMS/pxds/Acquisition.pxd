from libcpp cimport bool
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/Acquisition.h>" namespace "OpenMS":
    
    cdef cppclass Acquisition(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        Acquisition() nogil except +
        Acquisition(Acquisition) nogil except +

        bool operator==(Acquisition &rhs) nogil except +
        bool operator!=(Acquisition &rhs) nogil except +
        String getIdentifier() nogil except +
        void setIdentifier(const String &identifier) nogil except +

