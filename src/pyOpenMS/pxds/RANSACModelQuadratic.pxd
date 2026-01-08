from Types cimport *

cdef extern from "<OpenMS/ML/RANSAC/RANSACModelQuadratic.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelQuadratic:
       RansacModelQuadratic() except + nogil  # compiler
       RansacModelQuadratic(RansacModelQuadratic &) except + nogil  # compiler
