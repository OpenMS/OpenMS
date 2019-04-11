from Types cimport *
from String cimport *
from ProgressLogger cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>" namespace "OpenMS":

    cdef cppclass MRMDecoy(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        MRMDecoy()                       nogil except +
        MRMDecoy(MRMDecoy)               nogil except + # wrap-ignore

        void generateDecoys(TargetedExperiment& exp,
                            TargetedExperiment& dec,
                            String method,
                            double aim_decoy_fraction,
                            bool switchKR,
                            String decoy_tag,
                            int max_attempts,
                            double identity_threshold,
                            double precursor_mz_shift,
                            double product_mz_shift,
                            double product_mz_threshold,
                            libcpp_vector[String] fragment_types,
                            libcpp_vector[size_t] fragment_charges,
                            bool enable_specific_losses,
                            bool enable_unspecific_losses,
                            int round_decPow) nogil except +

        libcpp_vector[size_t] findFixedResidues(const String & sequence,
                                                bool keepN,
                                                bool keepC,
                                                const String & keep_const_pattern) nogil except +

