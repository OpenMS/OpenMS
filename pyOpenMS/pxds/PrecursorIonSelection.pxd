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
        PrecursorIonSelection() nogil except +
        PrecursorIonSelection(PrecursorIonSelection) nogil except +
        DoubleReal  getMaxScore() nogil except +
        void setMaxScore(DoubleReal & max_score) nogil except +
        void sortByTotalScore(FeatureMap[Feature] & features) nogil except +
        void getNextPrecursors(FeatureMap[Feature] & features, FeatureMap[Feature] & next_features, UInt number) nogil except +
        # TODO immutable types by reference
        # void getNextPrecursorsSeq(FeatureMap[Feature] & features, FeatureMap[Feature] & next_features, UInt number, DoubleReal & rt) nogil except +
        void getNextPrecursors(libcpp_vector[ int ] & solution_indices,
                               libcpp_vector[ IndexTriple ] & variable_indices,
                               libcpp_set[ int ] & measured_variables,
                               FeatureMap[Feature] & features,
                               FeatureMap[Feature] & new_features,
                               UInt step_size,
                               PSLPFormulation & ilp) nogil except +
        void rescore(FeatureMap[Feature] & features, libcpp_vector[ PeptideIdentification ] & new_pep_ids, libcpp_vector[ ProteinIdentification ] & prot_ids, PrecursorIonSelectionPreprocessing & preprocessed_db, bool check_meta_values) nogil except +
        void simulateRun(FeatureMap[Feature] & features, libcpp_vector[ PeptideIdentification ] & pep_ids, libcpp_vector[ ProteinIdentification ] & prot_ids, PrecursorIonSelectionPreprocessing & preprocessed_db, String path, MSExperiment[Peak1D, ChromatogramPeak] & experiment, String precursor_path) nogil except +
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
