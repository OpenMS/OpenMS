from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/MATH/MISC/CubicSpline2d.h>" namespace "OpenMS":
    
    cdef cppclass CubicSpline2d "OpenMS::CubicSpline2d":

        CubicSpline2d(libcpp_vector[ double ] x, libcpp_vector[ double ] y) except + nogil 
        CubicSpline2d(CubicSpline2d &) except + nogil  # compiler
        CubicSpline2d(libcpp_map[ double, double ] m) except + nogil 

        double eval(double x) except + nogil  # wrap-doc:Evaluates the cubic spline
        double derivatives(double x, unsigned order) except + nogil  # wrap-doc:Returns first, second or third derivative of cubic spline

