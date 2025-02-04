from Types cimport *

cdef extern from "<OpenMS/ML/RANSAC/RANSACModelLinear.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelLinear:
       RansacModelLinear() except + nogil  # compiler
       RansacModelLinear(RansacModelLinear &) except + nogil  # compiler
       # Other functions use iterators -> cannot be wrapped directly
