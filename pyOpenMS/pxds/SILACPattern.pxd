from Types cimport *
from SILACPoint cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACPattern.h>" namespace "OpenMS":
    
    cdef cppclass SILACPattern(SILACPoint) :
        # wrap-inherits:
        #  SILACPoint
        SILACPattern() nogil except +
        SILACPattern(SILACPattern) nogil except + #wrap-ignore
        # libcpp_vector[ SILACPoint ] points

