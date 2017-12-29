from Types cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/MRMMapping.h>" namespace "OpenMS":
    
    cdef cppclass MRMMapping(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MRMMapping() nogil except +
        MRMMapping(MRMMapping) nogil except + #wrap-ignore

        void mapExperiment(MSExperiment input_chromatograms, TargetedExperiment targeted_exp, MSExperiment& output) nogil except +

