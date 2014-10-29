from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/MATH/MISC/CubicSpline2d.h>" namespace "OpenMS":
    
    cdef cppclass CubicSpline2d "OpenMS::CubicSpline2d":

        CubicSpline2d(CubicSpline2d) nogil except + #wrap-ignore
        CubicSpline2d(libcpp_vector[ double ] x, libcpp_vector[ double ] y) nogil except +
        CubicSpline2d(libcpp_map[ double, double ] m) nogil except +

        double eval(double x) nogil except +
        double derivatives(double x, unsigned order) nogil except +

