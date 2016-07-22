from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from TargetedExperiment cimport *

# typedef boost::unordered_map<String, double> IonSeries;

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h>" namespace "OpenMS":
    
    cdef cppclass MRMIonSeries "OpenMS::MRMIonSeries":
        MRMIonSeries() nogil except +
        MRMIonSeries(MRMIonSeries) nogil except + #wrap-ignore

        ## Typedef boost::unordered_map
        ## libcpp_pair[ String, double ] getIon(IonSeries ionseries, String ionid) nogil except +
        ## libcpp_pair[ String, double ] annotateIon(IonSeries ionseries, double ProductMZ, double mz_threshold) nogil except +

        void annotateTransitionCV(ReactionMonitoringTransition & tr, String annotation) nogil except +

        void annotateTransition(ReactionMonitoringTransition & tr,
                                Peptide peptide, double
                                precursor_mz_threshold, double
                                product_mz_threshold, bool enable_reannotation,
                                libcpp_vector[ String ] fragment_types,
                                libcpp_vector[ size_t ] fragment_charges, 
                                bool enable_specific_losses, 
                                bool enable_unspecific_losses, int round_decPow) nogil except +

        ## IonSeries getIonSeries(AASequence sequence, size_t precursor_charge,
        ##                        libcpp_vector[ String ] fragment_types,
        ##                        libcpp_vector[ size_t ] fragment_charges, 
        ##                        bool enable_specific_losses,
        ##                        bool enable_unspecific_losses, int round_decPow) nogil except +

