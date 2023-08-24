from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSACModelQuadratic.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelQuadratic:
       RansacModelQuadratic() except + nogil  # compiler
       RansacModelQuadratic(RansacModelQuadratic &) except + nogil  # compiler
