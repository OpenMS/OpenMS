from ProgressLogger cimport *
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from Types cimport *

from TargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>" namespace "OpenMS":

    cdef cppclass MRMDecoy(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        MRMDecoy()                       nogil except +
        MRMDecoy(MRMDecoy)               nogil except + # wrap-ignore

        void generateDecoys(TargetedExperiment & exp,
                            TargetedExperiment & dec, String method, String decoy_tag,
                            double identity_threshold, int max_attempts, double mz_threshold, 
                            bool theoretical, double mz_shift, bool exclude_similar, 
                            double similarity_threshold, bool remove_CNterm_mods) nogil except +

        void restrictTransitions(TargetedExperiment & exp, int min_transitions, int max_transitions) nogil except +
