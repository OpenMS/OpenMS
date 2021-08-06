from Types cimport *
from DefaultParamHandler cimport *
from FeatureMap cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>" namespace "OpenMS":
    
    cdef cppclass ElutionModelFitter(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        ElutionModelFitter() nogil except + # wrap-doc:Helper class for fitting elution models to features
        ElutionModelFitter(ElutionModelFitter &) nogil except + # compiler
        void fitElutionModels(FeatureMap & features) nogil except + # wrap-doc:Fit models of elution profiles to all features (and validate them)

