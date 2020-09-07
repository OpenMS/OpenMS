from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSACModelLinear.h>" namespace "OpenMS::Math":

    cdef cppclass RansacModelLinear:
       RansacModelLinear() nogil except +
       RansacModelLinear(RansacModelLinear &) nogil except + # wrap-ignore
       # Other functions use iterators -> cannot be wrapped directly

