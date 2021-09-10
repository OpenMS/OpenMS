from Types cimport *
# from BaseLabeler cimport *
from DefaultParamHandler cimport *
from Param cimport *

cdef extern from "<OpenMS/SIMULATION/LABELING/SILACLabeler.h>" namespace "OpenMS":
    
    cdef cppclass SILACLabeler :
        # TODO BaseLabeler inheritance
        #  BaseLabeler
        SILACLabeler() nogil except + # wrap-doc:Simulate SILAC experiments
        SILACLabeler(SILACLabeler &) nogil except + # compiler

        void preCheck(Param & param) nogil except +
        ## void setUpHook(SimTypes::FeatureMapSimVector & ) nogil except +
        ## void postDigestHook(SimTypes::FeatureMapSimVector & ) nogil except +
        ## void postRTHook(SimTypes::FeatureMapSimVector & ) nogil except +
        ## void postDetectabilityHook(SimTypes::FeatureMapSimVector & ) nogil except +
        ## void postIonizationHook(SimTypes::FeatureMapSimVector & ) nogil except +
        ## void postRawMSHook(SimTypes::FeatureMapSimVector & ) nogil except +
        ## void postRawTandemMSHook(SimTypes::FeatureMapSimVector & , SimTypes::MSSimExperiment & ) nogil except +
        ## # POINTER # BaseLabeler * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Name of the model (needed by Factory)
