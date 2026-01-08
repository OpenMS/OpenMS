from Types cimport *
from DefaultParamHandler cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/FEATUREFINDER/ElutionModelFitter.h>" namespace "OpenMS":
    
    cdef cppclass ElutionModelFitter(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        ElutionModelFitter() except + nogil  # wrap-doc:Helper class for fitting elution models to features
        ElutionModelFitter(ElutionModelFitter &) except + nogil  # compiler
        void fitElutionModels(FeatureMap & features) except + nogil  # wrap-doc:Fit models of elution profiles to all features (and validate them)

