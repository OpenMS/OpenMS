from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from PeptideIdentification cimport *
from PSLPFormulation cimport *
from PrecursorIonSelectionPreprocessing cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>" namespace "OpenMS":
    
    cdef cppclass PrecursorIonSelection(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        # This class implements different precursor ion selection strategies
        PrecursorIonSelection() nogil except +
        PrecursorIonSelection(PrecursorIonSelection &) nogil except +

        double  getMaxScore() nogil except +
        void setMaxScore(double & max_score) nogil except +
        void sortByTotalScore(FeatureMap & features) nogil except + # wrap-doc:Sort features by total score
        void getNextPrecursors(FeatureMap & features, FeatureMap & next_features, UInt number) nogil except +
            # wrap-doc:
                #   Returns features with highest score for MS/MS
                #   -----
                #   :param features: FeatureMap with all possible precursors
                #   :param next_features: FeatureMap with next precursors
                #   :param number: Number of features to be reported

        # TODO immutable types by reference
        # void getNextPrecursorsSeq(FeatureMap & features, FeatureMap & next_features, UInt number, double & rt) nogil except +

        void getNextPrecursors(libcpp_vector[ int ] & solution_indices,
                               libcpp_vector[ IndexTriple ] & variable_indices,
                               libcpp_set[ int ] & measured_variables,
                               FeatureMap & features,
                               FeatureMap & new_features,
                               UInt step_size,
                               PSLPFormulation & ilp) nogil except +
        void rescore(FeatureMap & features, libcpp_vector[ PeptideIdentification ] & new_pep_ids, libcpp_vector[ ProteinIdentification ] & prot_ids, PrecursorIonSelectionPreprocessing & preprocessed_db, bool check_meta_values) nogil except +
            # wrap-doc:
                #   Change scoring of features using peptide identifications from all spectra
                #   -----
                #   :param features: FeatureMap with all possible precursors
                #   :param new_pep_ids: Peptide identifications
                #   :param prot_ids: Protein identifications
                #   :param preprocessed_db: Information from preprocessed database
                #   :param check_meta_values: True if the FeatureMap should be checked for the presence of required meta values

        void simulateRun(FeatureMap & features, libcpp_vector[ PeptideIdentification ] & pep_ids, libcpp_vector[ ProteinIdentification ] & prot_ids, PrecursorIonSelectionPreprocessing & preprocessed_db, String path, MSExperiment & experiment, String precursor_path) nogil except +
            # wrap-doc:
                #   Simulate the iterative precursor ion selection
                #   -----
                #   :param features: FeatureMap with all possible precursors
                #   :param new_pep_ids: Peptide identifications
                #   :param prot_ids: Protein identifications
                #   :param preprocessed_db: Information from preprocessed database
                #   :param step_size: Number of MS/MS spectra considered per iteration
                #   :param path: Path to output file

        void setLPSolver(SOLVER solver) nogil except +
        SOLVER getLPSolver() nogil except +
        void reset() nogil except +
        # TODO not implemented but part of the API
        # libcpp_map[ String, libcpp_set[ String ] ]  getPeptideProteinCounter() nogil except +

cdef extern from "<OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>" namespace "OpenMS::PrecursorIonSelection":
    cdef enum PrecursorIonSelection_Type "OpenMS::PrecursorIonSelection::Type":
        #wrap-attach:
        #    PrecursorIonSelection
        IPS
        ILP_IPS
        SPS
        UPSHIFT
        DOWNSHIFT
        DEX

# cdef extern from "<OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>" namespace "OpenMS::PrecursorIonSelection":
#     
#     cdef cppclass SeqTotalScoreMore :
#         SeqTotalScoreMore(SeqTotalScoreMore) nogil except + #wrap-ignore
#         bool operator()(Feature & left, Feature & right) nogil except + # wrap-cast:evaluate
# 
# cdef extern from "<OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>" namespace "OpenMS::PrecursorIonSelection":
#     
#     cdef cppclass TotalScoreMore :
#         TotalScoreMore(TotalScoreMore) nogil except + #wrap-ignore
#         bool operator()(Feature & left, Feature & right) nogil except + # wrap-cast:evaluate
# 
