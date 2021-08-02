from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from TargetedExperiment cimport *

# typedef boost::unordered_map<String, double> IonSeries;

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h>" namespace "OpenMS":
    
    cdef cppclass MRMIonSeries "OpenMS::MRMIonSeries":
        MRMIonSeries() nogil except +
        MRMIonSeries(MRMIonSeries &) nogil except + # compiler

        ## Typedef boost::unordered_map
        ## libcpp_pair[ String, double ] getIon(IonSeries ionseries, String ionid) nogil except +
        ## libcpp_pair[ String, double ] annotateIon(IonSeries ionseries, double ProductMZ, double mz_threshold) nogil except +

        void annotateTransitionCV(ReactionMonitoringTransition & tr, String annotation) nogil except +
            # wrap-doc:
                #   Annotates transition with CV terms
                #   -----
                #   :param tr: The transition to annotate
                #   :param annotation: The fragment ion annotation

        void annotateTransition(ReactionMonitoringTransition & tr,
                                Peptide peptide, double
                                precursor_mz_threshold, double
                                product_mz_threshold, bool enable_reannotation,
                                libcpp_vector[ String ] fragment_types,
                                libcpp_vector[ size_t ] fragment_charges, 
                                bool enable_specific_losses, 
                                bool enable_unspecific_losses, int round_decPow) nogil except +
            # wrap-doc:
                #   Annotates transition
                #   -----
                #   :param tr: The transition to annotate
                #   :param peptide: The corresponding peptide
                #   :param precursor_mz_threshold: The m/z threshold for annotation of the precursor ion
                #   :param product_mz_threshold: The m/z threshold for annotation of the fragment ion
                #   :param enable_reannotation: Whether the original (e.g. SpectraST) annotation should be used or reannotation should be conducted
                #   :param fragment_types: The fragment ion types for reannotation
                #   :param fragment_charges: The fragment ion charges for reannotation
                #   :param enable_specific_losses: Whether specific neutral losses should be considered
                #   :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
                #   :param round_decPow: Round precursor and product m/z values to decimal power (default: -4)

        ## IonSeries getIonSeries(AASequence sequence, size_t precursor_charge,
        ##                        libcpp_vector[ String ] fragment_types,
        ##                        libcpp_vector[ size_t ] fragment_charges, 
        ##                        bool enable_specific_losses,
        ##                        bool enable_unspecific_losses, int round_decPow) nogil except +

