from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSACModelQuadratic.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelQuadratic:
       RansacModelQuadratic() nogil except + # compiler
       RansacModelQuadratic(RansacModelQuadratic &) nogil except + # compiler
