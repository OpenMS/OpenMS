from libc.string cimport const_char


cdef extern from "<OpenMS/DATASTRUCTURES/String.h>" namespace "OpenMS":

    cdef cppclass String:
         String() nogil except +
         String(String) nogil except +  # wrap-ignore
         String(char *) nogil except +
         const_char * c_str() nogil except +
