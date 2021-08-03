from LevMarqFitter1D cimport *
from InterpolationModel cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass EmgFitter1D(LevMarqFitter1D):
        # wrap-inherits:
        #  LevMarqFitter1D
        EmgFitter1D() nogil except + # wrap-doc:Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization
        EmgFitter1D(EmgFitter1D &) nogil except +
        # float fit1d(libcpp_vector[Peak1D] range_, InterpolationModel * & model) nogil except + # wrap-ignore
        # Fitter1D * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Name of the model (needed by Factory)

