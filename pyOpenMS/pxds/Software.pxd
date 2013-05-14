from String cimport *
from Types cimport *

cdef extern from "<OpenMS/METADATA/Software.h>" namespace "OpenMS":

    cdef cppclass Software:

        Software()   nogil except +
        Software(Software) nogil except + # wrap-ignore

        String getName() nogil except +
        String getVersion() nogil except +

        void setName(String) nogil except +
        void setVersion(String) nogil except +



