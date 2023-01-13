from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>" namespace "OpenMS":

    cdef cppclass FalseDiscoveryRate(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        FalseDiscoveryRate() nogil except +
        # private
        FalseDiscoveryRate(FalseDiscoveryRate &) nogil except + #wrap-ignore

        void apply(libcpp_vector[PeptideIdentification] & forward_ids, libcpp_vector[PeptideIdentification] & reverse_ids) nogil except + 
        void apply(libcpp_vector[PeptideIdentification] & id) nogil except + 
        void apply(libcpp_vector[ProteinIdentification] & forward_ids, libcpp_vector[ProteinIdentification] & reverse_ids) nogil except + 
        void apply(libcpp_vector[ProteinIdentification] & id) nogil except + 

        void applyEstimated(libcpp_vector[ProteinIdentification] & ids) nogil except +
        double applyEvaluateProteinIDs(libcpp_vector[ProteinIdentification] & ids, double pepCutoff, UInt fpCutoff, double diffWeight) nogil except +
        double applyEvaluateProteinIDs(ProteinIdentification& ids, double pepCutoff, UInt fpCutoff, double diffWeight) nogil except +

        # simpler reimplementation of the apply function above for PSMs. With charge and identifier info from run_info
        void applyBasic(libcpp_vector[ProteinIdentification] & run_info, libcpp_vector[PeptideIdentification] & ids) nogil except +
        # simpler reimplementation of the apply function above for PSMs or peptides. Use charge=0 and identifier="" for "all".
        void applyBasic(libcpp_vector[PeptideIdentification] & ids, bool higher_score_better, int charge, String identifier, bool only_best_per_pep);
    
        # simpler reimplementation of the apply function above for peptides in ConsensusMaps.
        void applyBasic(ConsensusMap & cmap, bool use_unassigned_peptides) nogil except +
        # simpler reimplementation of the apply function above for proteins.
        void applyBasic(ProteinIdentification & id, bool groups_too) nogil except +
        # applies a picked protein FDR (TODO explain/ref)
        void applyPickedProteinFDR(ProteinIdentification & id, String & decoy_string, bool decoy_prefix, bool groups_too) nogil except +

        # calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
        # if fp_cutoff = 0, it will calculate the full AUC
        double rocN(libcpp_vector[PeptideIdentification] & ids, Size fp_cutoff) nogil except +
        # calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
        # if fp_cutoff = 0, it will calculate the full AUC
        double rocN(ConsensusMap& ids, Size fp_cutoff, bool include_unassigned_peptides) nogil except +
        # calculates the AUC until the first fp_cutoff False positive pep IDs for a specific run with identifier
        # if fp_cutoff = 0, it will calculate the full AUC
        double rocN(ConsensusMap& ids, Size fp_cutoff, const String& identifier, bool include_unassigned_peptides) nogil except +

