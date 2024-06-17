from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from EmgFitter1D cimport *
from EmgModel cimport *
from GaussFilter cimport *
from MRMFeature cimport *
from MRMTransitionGroup cimport *
from DPosition cimport *

cdef extern from "<OpenMS/FEATUREFINDER/EmgScoring.h>" namespace "OpenMS":
    
    cdef cppclass EmgScoring "OpenMS::EmgScoring":
        EmgScoring() except + nogil  # wrap-doc:Helps in scoring of an elution peak using an exponentially modified gaussian distribution model
        EmgScoring(EmgScoring &) except + nogil  # compiler

        void setFitterParam(Param param) except + nogil  # TODO
        Param getDefaults() except + nogil  # TODO
        # TEMPLATE # double calcElutionFitScore(MRMFeature &mrmfeature, MRMTransitionGroup[ SpectrumType, TransitionT ] &transition_group) except + nogil 
        double elutionModelFit(libcpp_vector[DPosition2] current_section, bool smooth_data) except + nogil 

