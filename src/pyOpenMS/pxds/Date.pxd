from String cimport *
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Date.h>" namespace "OpenMS":

    cdef cppclass Date:

        Date() nogil except +
        Date(Date) nogil except + # wrap-ignore

        void set(String & date) nogil except +
        # void set(UInt month, UInt day, UInt year);

        Date today() nogil except +
        String get()  nogil except +
        # void get(UInt & month, UInt & day, UInt & year) nogil except +

        # Sets the undefined date: 00/00/0000
        void clear() nogil except +
