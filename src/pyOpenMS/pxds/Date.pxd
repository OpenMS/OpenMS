from String cimport *
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Date.h>" namespace "OpenMS":

    cdef cppclass Date:

        Date() except + nogil 
        Date(Date &) except + nogil 

        void set(const String & date) except + nogil 
        # void set(UInt month, UInt day, UInt year);

        Date today() except + nogil 
        String get()  except + nogil 
        # void get(UInt & month, UInt & day, UInt & year) except + nogil 

        # Sets the undefined date: 00/00/0000
        void clear() except + nogil 
