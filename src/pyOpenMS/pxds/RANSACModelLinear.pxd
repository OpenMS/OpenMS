from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSACModelLinear.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelLinear:
       RansacModelLinear() nogil except + # compiler
       RansacModelLinear(RansacModelLinear &) nogil except + # compiler
       # Other functions use iterators -> cannot be wrapped directly
