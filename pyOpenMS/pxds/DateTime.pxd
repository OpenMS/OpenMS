from String cimport *
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DateTime.h>" namespace "OpenMS":

    cdef cppclass DateTime:

        DateTime()   nogil except +
        DateTime(DateTime) nogil except + # wrap-ignore

        String getDate() nogil except +
        String getTime() nogil except +


cdef extern from "<OpenMS/DATASTRUCTURES/DateTime.h>" namespace "OpenMS::DateTime":

    DateTime now() # wrap-attach:DateTime




