from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from EmgFitter1D cimport *
from EmgModel cimport *
from GaussFilter cimport *
from MRMFeature cimport *
from MRMTransitionGroup cimport *
from DPosition cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgScoring.h>" namespace "OpenMS":
    
    cdef cppclass EmgScoring "OpenMS::EmgScoring":
        EmgScoring() nogil except + # wrap-doc:Helps in scoring of an elution peak using an exponentially modified gaussian distribution model
        EmgScoring(EmgScoring &) nogil except + # compiler

        void setFitterParam(Param param) nogil except + # TODO
        Param getDefaults() nogil except + # TODO
        # TEMPLATE # double calcElutionFitScore(MRMFeature &mrmfeature, MRMTransitionGroup[ SpectrumType, TransitionT ] &transition_group) nogil except +
        double elutionModelFit(libcpp_vector[DPosition2] current_section, bool smooth_data) nogil except +

