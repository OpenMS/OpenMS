from Types cimport *
from DefaultParamHandler cimport *
from FeatureMap cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>" namespace "OpenMS":
    
    cdef cppclass ElutionModelFitter(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        ElutionModelFitter() nogil except +
        ElutionModelFitter(ElutionModelFitter) nogil except + #wrap-ignore
        void fitElutionModels(FeatureMap & features) nogil except +

