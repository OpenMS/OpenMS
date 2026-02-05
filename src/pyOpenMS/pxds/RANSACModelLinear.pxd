from Types cimport *

cdef extern from "<OpenMS/ML/RANSAC/RANSACModelLinear.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelLinear:
        # wrap-doc:
        #  Implementation of a linear RANSAC model fit
       RansacModelLinear() except + nogil  # compiler
       RansacModelLinear(RansacModelLinear &) except + nogil  # compiler
       # Other functions use iterators -> cannot be wrapped directly
