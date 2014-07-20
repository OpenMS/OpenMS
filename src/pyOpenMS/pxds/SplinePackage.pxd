from Types cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>" namespace "OpenMS":
    
    cdef cppclass SplinePackage "OpenMS::SplinePackage":

        SplinePackage(libcpp_vector[double] mz, libcpp_vector[double] intensity, double scaling) nogil except +
        SplinePackage(SplinePackage) nogil except + #wrap-ignore

        double getMzMin() nogil except +
        double getMzMax() nogil except +
        double getMzStepWidth() nogil except +
        bool isInPackage(double mz) nogil except +
        double eval(double mz) nogil except +

