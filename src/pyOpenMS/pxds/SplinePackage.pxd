from Types cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>" namespace "OpenMS":
    
    cdef cppclass SplinePackage "OpenMS::SplinePackage":

        SplinePackage(libcpp_vector[double] pos, libcpp_vector[double] intensity) nogil except +
        SplinePackage(SplinePackage) nogil except + #wrap-ignore

        double getPosMin() nogil except +
        double getPosMax() nogil except +
        double getPosStepWidth() nogil except +
        bool isInPackage(double pos) nogil except +
        double eval(double pos) nogil except +

