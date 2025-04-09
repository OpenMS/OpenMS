from String cimport *
from Types cimport *

cdef extern from "<OpenMS/METADATA/Software.h>" namespace "OpenMS":

    cdef cppclass Software:

        Software() except + nogil 
        Software(Software &) except + nogil 

        String getName() except + nogil  # wrap-doc:Returns the name of the software
        String getVersion() except + nogil  # wrap-doc:Returns the software version

        void setName(String) except + nogil  # wrap-doc:Sets the name of the software
        void setVersion(String) except + nogil  # wrap-doc:Sets the software version
