from Types cimport *
from TargetedExperiment cimport *
from ProgressLogger cimport *
from MRMIonSeries cimport *
from ModificationsDB cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>" namespace "OpenMS":
    
    cdef cppclass MRMAssay(ProgressLogger) :
        # wrap-inherits:
        #  ProgressLogger

        MRMAssay() nogil except +
        MRMAssay(MRMAssay) nogil except + #wrap-ignore

        void reannotateTransitions(TargetedExperiment & exp, 
                                   double precursor_mz_threshold,
                                   double product_mz_threshold,
                                   libcpp_vector[ String ] fragment_types, 
                                   libcpp_vector[ size_t ] fragment_charges, 
                                   bool enable_specific_losses, 
                                   bool enable_unspecific_losses,
                                   int round_decPow) nogil except +

        void restrictTransitions(TargetedExperiment & exp, 
                                 double lower_mz_limit,
                                 double upper_mz_limit,
                                 libcpp_vector[ libcpp_pair[ double, double ] ] swathes) nogil except +

        void detectingTransitions(TargetedExperiment & exp,
                                  int min_transitions, int max_transitions) nogil except +

        void uisTransitions(TargetedExperiment & exp, 
                            libcpp_vector[ String ] fragment_types,
                            libcpp_vector[ size_t ] fragment_charges,
                            bool enable_specific_losses,
                            bool enable_unspecific_losses,
                            bool enable_ms2_precursors,
                            double mz_threshold,
                            libcpp_vector[ libcpp_pair[ double, double ] ] swathes,
                            int round_decPow,
                            size_t max_num_alternative_localizations,
                            int shuffle_seed) nogil except +

