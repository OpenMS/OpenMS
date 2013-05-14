from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/Element.h>" namespace "OpenMS":

    cdef cppclass Element:

        Element() nogil except +
        Element(Element) nogil except + # wrap-ignore

