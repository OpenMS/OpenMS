from Fitter1D cimport *

cdef extern from "<OpenMS/FEATUREFINDER/LevMarqFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass LevMarqFitter1D(Fitter1D):
        # wrap-ignore
        # no-pxd-import
        # wrap-doc:
        # Abstract class for 1D-model fitter using Levenberg-Marquardt algorithm for parameter optimization
        
        LevMarqFitter1D() except + nogil 
        LevMarqFitter1D(LevMarqFitter1D &) except + nogil 
