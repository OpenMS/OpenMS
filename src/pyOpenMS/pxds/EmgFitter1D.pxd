from LevMarqFitter1D cimport *
from InterpolationModel cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FEATUREFINDER/EmgFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass EmgFitter1D(LevMarqFitter1D):
        # wrap-inherits:
        #  LevMarqFitter1D
        EmgFitter1D() except + nogil  # wrap-doc:Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization
        EmgFitter1D(EmgFitter1D &) except + nogil 
        # float fit1d(libcpp_vector[Peak1D] range_, InterpolationModel * & model) except + nogil  # wrap-ignore
        # Fitter1D * create() except + nogil 
      
