#from libcpp.string cimport *
#from libcpp.vector cimport *
#from DataValue cimport *
#from String cimport *
from libc.string cimport const_char


cdef extern from "<OpenMS/DATASTRUCTURES/String.h>" namespace "OpenMS":

    cdef cppclass String:
        # we have converters for this, do not wrap the class itself !
        # wrap-ignore
         String() nogil except +
         String(String) nogil except +
         String(char *) nogil except +
         const_char * c_str() nogil except +
