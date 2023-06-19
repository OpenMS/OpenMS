from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/BSpline2d.h>" namespace "OpenMS":
    
    cdef cppclass BSpline2d:

        BSpline2d(libcpp_vector[double] x, libcpp_vector[double] y, 
                  double wave_length, BoundaryCondition boundary_condition, 
                  Size num_nodes) nogil except +

        bool solve(libcpp_vector[double] y) nogil except + # wrap-doc:Solve the spline curve for a new set of y values. Returns false if the solution fails

        double eval(double x) nogil except + # wrap-doc:Returns the evaluation of the smoothed curve at a particular x value. If current state is not ok(), returns zero

        double derivative(double x) nogil except + # wrap-doc:Returns the first derivative of the spline curve at the given position x. Returns zero if the current state is not ok()

        bool ok() nogil except + # wrap-doc:Returns whether the spline fit was successful

        void debug(bool enable) nogil except + # wrap-doc:Enable or disable debug messages from the B-spline library

cdef extern from "<OpenMS/MATH/MISC/BSpline2d.h>" namespace "OpenMS::BSpline2d":

    cdef enum BoundaryCondition:
        BC_ZERO_ENDPOINTS, BC_ZERO_FIRST, BC_ZERO_SECOND
