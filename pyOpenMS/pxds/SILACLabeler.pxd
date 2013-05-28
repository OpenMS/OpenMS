from Types cimport *
# from BaseLabeler cimport *
from DefaultParamHandler cimport *
from Param cimport *

cdef extern from "<OpenMS/SIMULATION/LABELING/SILACLabeler.h>" namespace "OpenMS":
    
    cdef cppclass SILACLabeler :
        # TODO BaseLabeler inheritance
        #  BaseLabeler
        SILACLabeler() nogil except +
        SILACLabeler(SILACLabeler) nogil except + #wrap-ignore
        void preCheck(Param & param) nogil except +
        ## void setUpHook(FeatureMapSimVector & ) nogil except +
        ## void postDigestHook(FeatureMapSimVector & ) nogil except +
        ## void postRTHook(FeatureMapSimVector & ) nogil except +
        ## void postDetectabilityHook(FeatureMapSimVector & ) nogil except +
        ## void postIonizationHook(FeatureMapSimVector & ) nogil except +
        ## void postRawMSHook(FeatureMapSimVector & ) nogil except +
        ## void postRawTandemMSHook(FeatureMapSimVector & , MSSimExperiment & ) nogil except +
        ## # POINTER # BaseLabeler * create() nogil except +
        String getProductName() nogil except +

