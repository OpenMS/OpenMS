from libcpp.vector cimport vector as libcpp_vector
from IsotopeCluster cimport *
from DefaultParamHandler cimport *
from Feature cimport *
# from BasicStatistics cimport *
# from FeatureFinderDefs cimport *
from Peak1D cimport *

ctypedef libcpp_vector[Peak1D] RawDataArrayType

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>" namespace "OpenMS":
    
    cdef cppclass Fitter1D(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        Fitter1D() nogil except + # wrap-doc:Abstract base class for all 1D-dimensional model fitter
        Fitter1D(Fitter1D &) nogil except +

        # QualityType fit1d(RawDataArrayType &, InterpolationModel *&) nogil except +
        void registerChildren() nogil except + # wrap-doc:Register all derived classes here
