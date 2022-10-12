from String cimport *
from Types cimport *

cdef extern from "<OpenMS/METADATA/Software.h>" namespace "OpenMS":

    cdef cppclass Software:

        Software() nogil except +
        Software(Software &) nogil except +

        String getName() nogil except + # wrap-doc:\n Returns the name of the software
        String getVersion() nogil except + # wrap-doc:\n Returns the software version

        void setName(String) nogil except + # wrap-doc:\n Sets the name of the software
        void setVersion(String) nogil except + # wrap-doc:\n Sets the software version
