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
                            double mz_shift, bool exclude_similar, 
                            double similarity_threshold, bool remove_CNterm_mods,
                            double precursor_mass_shift, libcpp_vector[String] fragment_types,
                            libcpp_vector[size_t] fragment_charges, bool enable_specific_losses,
                            bool enable_unspecific_losses, bool remove_unannotated,
                            int round_decPow) nogil except +

