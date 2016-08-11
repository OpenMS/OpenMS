from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSACModelQuadratic.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelQuadratic:
       RansacModelQuadratic() nogil except +
       RansacModelQuadratic(RansacModelQuadratic &) nogil except + # wrap-ignore

