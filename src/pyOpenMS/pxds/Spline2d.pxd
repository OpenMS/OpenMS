from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/Spline2d.h>" namespace "OpenMS":
    
    cdef cppclass Spline2d[ValType]:
        # wrap-instances:
        #   Spline2d := Spline2d[double]

        Spline2d() nogil except + #wrap-ignore
        Spline2d(Spline2d) nogil except + #wrap-ignore

        Spline2d(unsigned degree, libcpp_vector[ValType] x, libcpp_vector[double] y) nogil except +

        # ValType eval(ValType x) nogil except +
        # ValType derivatives(ValType x, unsigned order) nogil except +

