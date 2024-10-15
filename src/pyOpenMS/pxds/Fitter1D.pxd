from libcpp.vector cimport vector as libcpp_vector
from IsotopeCluster cimport *
from DefaultParamHandler cimport *
from Feature cimport *
# from BasicStatistics cimport *
# from FeatureFinderDefs cimport *
from Peak1D cimport *

ctypedef libcpp_vector[Peak1D] RawDataArrayType

cdef extern from "<OpenMS/FEATUREFINDER/Fitter1D.h>" namespace "OpenMS":
    
    cdef cppclass Fitter1D(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        Fitter1D() except + nogil  # wrap-doc:Abstract base class for all 1D-dimensional model fitter
        Fitter1D(Fitter1D &) except + nogil 

        # QualityType fit1d(RawDataArrayType &, InterpolationModel *&) except + nogil 
