from Types cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>" namespace "OpenMS":
    
    cdef cppclass SplinePackage "OpenMS::SplinePackage":

        SplinePackage(libcpp_vector[double] pos, libcpp_vector[double] intensity) nogil except +
        SplinePackage(SplinePackage &) nogil except + # compiler

        double getPosMin() nogil except + # wrap-doc:Returns the minimum position for which the spline fit is valid
        double getPosMax() nogil except + # wrap-doc:Returns the maximum position for which the spline fit is valid
        double getPosStepWidth() nogil except + # wrap-doc:Returns a sensible position step width for the package
        bool isInPackage(double pos) nogil except + # wrap-doc:Returns true if position in [posMin:posMax] interval else false
        double eval(double pos) nogil except + # wrap-doc:Returns interpolated intensity position `pos`
