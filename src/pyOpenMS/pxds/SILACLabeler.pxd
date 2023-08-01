from Types cimport *
# from BaseLabeler cimport *
from DefaultParamHandler cimport *
from Param cimport *

cdef extern from "<OpenMS/SIMULATION/LABELING/SILACLabeler.h>" namespace "OpenMS":
    
    cdef cppclass SILACLabeler :
        # TODO BaseLabeler inheritance
        #  BaseLabeler
        SILACLabeler() except + nogil  # wrap-doc:Simulate SILAC experiments
        SILACLabeler(SILACLabeler &) except + nogil  # compiler

        void preCheck(Param & param) except + nogil 
        ## void setUpHook(SimTypes::FeatureMapSimVector & ) except + nogil 
        ## void postDigestHook(SimTypes::FeatureMapSimVector & ) except + nogil 
        ## void postRTHook(SimTypes::FeatureMapSimVector & ) except + nogil 
        ## void postDetectabilityHook(SimTypes::FeatureMapSimVector & ) except + nogil 
        ## void postIonizationHook(SimTypes::FeatureMapSimVector & ) except + nogil 
        ## void postRawMSHook(SimTypes::FeatureMapSimVector & ) except + nogil 
        ## void postRawTandemMSHook(SimTypes::FeatureMapSimVector & , SimTypes::MSSimExperiment & ) except + nogil 
        ## # POINTER # BaseLabeler * create() except + nogil 
        String getProductName() except + nogil  # wrap-doc:Name of the model (needed by Factory)
