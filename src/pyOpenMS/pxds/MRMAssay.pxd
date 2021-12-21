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
        MRMAssay(MRMAssay &) nogil except + # compiler

        void reannotateTransitions(TargetedExperiment & exp, 
                                   double precursor_mz_threshold,
                                   double product_mz_threshold,
                                   libcpp_vector[ String ] fragment_types, 
                                   libcpp_vector[ size_t ] fragment_charges, 
                                   bool enable_specific_losses, 
                                   bool enable_unspecific_losses,
                                   int round_decPow) nogil except +
            # wrap-doc:
                #   Annotates and filters transitions in a TargetedExperiment
                #   -----
                #   :param exp: The input, unfiltered transitions
                #   :param precursor_mz_threshold: The precursor m/z threshold in Th for annotation
                #   :param product_mz_threshold: The product m/z threshold in Th for annotation
                #   :param fragment_types: The fragment types to consider for annotation
                #   :param fragment_charges: The fragment charges to consider for annotation
                #   :param enable_specific_losses: Whether specific neutral losses should be considered
                #   :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
                #   :param round_decPow: Round product m/z values to decimal power (default: -4)

        void restrictTransitions(TargetedExperiment & exp, 
                                 double lower_mz_limit,
                                 double upper_mz_limit,
                                 libcpp_vector[ libcpp_pair[ double, double ] ] swathes) nogil except +
            # wrap-doc:
                #   Restrict and filter transitions in a TargetedExperiment
                #   -----
                #   :param exp: The input, unfiltered transitions
                #   :param lower_mz_limit: The lower product m/z limit in Th
                #   :param upper_mz_limit: The upper product m/z limit in Th
                #   :param swathes: The swath window settings (to exclude fragment ions falling into the precursor isolation window)

        void detectingTransitions(TargetedExperiment & exp,
                                  int min_transitions,
                                  int max_transitions) nogil except +
            # wrap-doc:
                #   Select detecting fragment ions
                #   -----
                #   :param exp: The input, unfiltered transitions
                #   :param min_transitions: The minimum number of transitions required per assay
                #   :param max_transitions: The maximum number of transitions required per assay

        void filterMinMaxTransitionsCompound(TargetedExperiment & exp,
                                          int min_transitions,
                                          int max_transitions) nogil except +
            # wrap-doc:
                #   Filters target and decoy transitions by intensity, only keeping the top N transitions
                #   -----
                #   :param exp: The transition list which will be filtered
                #   :param min_transitions: The minimum number of transitions required per assay (targets only)
                #   :param max_transitions: The maximum number of transitions allowed per assay

        void filterUnreferencedDecoysCompound(TargetedExperiment & exp) nogil except +
            # wrap-doc:
                #   Filters decoy transitions, which do not have respective target transition
                #   based on the transitionID.
                #   -----
                #   References between targets and decoys will be constructed based on the transitionsID
                #   and the "_decoy_" string. For example:
                #   -----
                #   target: 84_CompoundName_[M+H]+_88_22
                #   decoy: 84_CompoundName_decoy_[M+H]+_88_22
                #   -----
                #   :param exp: The transition list which will be filtered
        
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
            # wrap-doc:
                #   Annotate UIS / site-specific transitions
                #   -----
                #   Performs the following actions:
                #   -----
                #   - Step 1: For each peptide, compute all theoretical alternative peptidoforms; see transitions generateTargetInSilicoMap_()
                #   - Step 2: Generate target identification transitions; see generateTargetAssays_()
                #   -----
                #   - Step 3a: Generate decoy sequences that share peptidoform properties with targets; see generateDecoySequences_()
                #   - Step 3b: Generate decoy in silico peptide map containing theoretical transition; see generateDecoyInSilicoMap_()
                #   - Step 4: Generate decoy identification transitions; see generateDecoyAssays_()
                #   -----
                #   The IPF algorithm uses the concept of "identification transitions" that
                #   are used to discriminate different peptidoforms, these are generated in
                #   this function.  In brief, the algorithm takes the existing set of
                #   peptides and transitions and then appends these "identification
                #   transitions" for targets and decoys. The novel transitions are set to be
                #   non-detecting and non-quantifying and are annotated with the set of
                #   peptidoforms to which they map.
                #   -----
                #   :param exp: The input, unfiltered transitions
                #   :param fragment_types: The fragment types to consider for annotation
                #   :param fragment_charges: The fragment charges to consider for annotation
                #   :param enable_specific_losses: Whether specific neutral losses should be considered
                #   :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
                #   :param enable_ms2_precursors: Whether MS2 precursors should be considered
                #   :param mz_threshold: The product m/z threshold in Th for annotation
                #   :param swathes: The swath window settings (to exclude fragment ions falling
                #   :param round_decPow: Round product m/z values to decimal power (default: -4)
                #   :param max_num_alternative_localizations: Maximum number of allowed peptide sequence permutations
                #   :param shuffle_seed: Set seed for shuffle (-1: select seed based on time)
                #   :param disable_decoy_transitions: Whether to disable generation of decoy UIS transitions

