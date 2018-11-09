from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/BSpline2d.h>" namespace "OpenMS":
    
    cdef cppclass BSpline2d:

        BSpline2d(libcpp_vector[double] x, libcpp_vector[double] y, 
                  double wave_length, BoundaryCondition boundary_condition, 
                  Size num_nodes) nogil except +

        bool solve(libcpp_vector[double] y) nogil except +

        double eval(double x) nogil except +

        double derivative(double x) nogil except +

        bool ok() nogil except +

        void debug(bool enable) nogil except +

        double derivatives(double x, unsigned order) nogil except +

cdef extern from "<OpenMS/MATH/MISC/BSpline2d.h>" namespace "OpenMS::BSpline2d":

    cdef enum BoundaryCondition:
        BC_ZERO_ENDPOINTS, BC_ZERO_FIRST, BC_ZERO_SECOND

