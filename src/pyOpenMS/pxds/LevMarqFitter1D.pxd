from Fitter1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/LevMarqFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass LevMarqFitter1D(Fitter1D):
        # wrap-ignore
        LevMarqFitter1D() nogil except +
        LevMarqFitter1D(LevMarqFitter1D) nogil except +

