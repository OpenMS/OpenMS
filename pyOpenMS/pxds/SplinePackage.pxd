from Types cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>" namespace "OpenMS":
    
    cdef cppclass SplinePackage "OpenMS::SplinePackage":
        SplinePackage() nogil except +
        SplinePackage(SplinePackage) nogil except + #wrap-ignore
        double mzMin_
        double mzMax_
	double mzStepWidth_
	Spline2d<double> spline_

