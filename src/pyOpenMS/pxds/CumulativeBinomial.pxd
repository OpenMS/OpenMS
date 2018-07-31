from Types cimport *

cdef extern from "<OpenMS/MATH/STATISTICS/CumulativeBinomial.h>" namespace "OpenMS::Math":
    
    cdef cppclass CumulativeBinomial "OpenMS::Math::CumulativeBinomial":
        CumulativeBinomial(CumulativeBinomial) nogil except + #wrap-ignore

        double compute(Size n, Size k, double p) nogil except +

