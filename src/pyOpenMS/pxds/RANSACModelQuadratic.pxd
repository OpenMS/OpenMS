from Types cimport *

cdef extern from "<OpenMS/ML/RANSAC/RANSACModelQuadratic.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelQuadratic:
        # wrap-doc:
        #  Implementation of a quadratic RANSAC model fit
       RansacModelQuadratic() except + nogil  # compiler
       RansacModelQuadratic(RansacModelQuadratic &) except + nogil  # compiler
