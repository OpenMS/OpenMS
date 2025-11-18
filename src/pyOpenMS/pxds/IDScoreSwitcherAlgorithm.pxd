from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from DefaultParamHandler cimport *
from String cimport *
from ConsensusMap cimport *
from PeptideIdentificationList cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>" namespace "OpenMS":

    cdef cppclass IDScoreSwitcherAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        IDScoreSwitcherAlgorithm() except + nogil
        IDScoreSwitcherAlgorithm(IDScoreSwitcherAlgorithm&) except + nogil

        bool isScoreType(const String& score_name, IDScoreSwitcherAlgorithm_ScoreType score_type) except + nogil
        @staticmethod
        IDScoreSwitcherAlgorithm_ScoreType toScoreTypeEnum(String score_type) except + nogil
        bool isScoreTypeHigherBetter(IDScoreSwitcherAlgorithm_ScoreType score_type) except + nogil
        libcpp_vector[String] getScoreNames() except + nogil
        void switchToGeneralScoreType(PeptideIdentificationList& pep_ids,
                                      IDScoreSwitcherAlgorithm_ScoreType score_type,
                                      Size& counter) except + nogil
        void switchToGeneralScoreType(ConsensusMap& cmap,
                                      IDScoreSwitcherAlgorithm_ScoreType score_type,
                                      Size& counter,
                                      bool unassigned_peptides_too) except + nogil
        void determineScoreNameOrientationAndType(const PeptideIdentificationList& pep_ids,
                                                  String& name,
                                                  bool& higher_better,
                                                  IDScoreSwitcherAlgorithm_ScoreType& score_type) except + nogil
        void determineScoreNameOrientationAndType(const ConsensusMap& cmap,
                                                  String& name,
                                                  bool& higher_better,
                                                  IDScoreSwitcherAlgorithm_ScoreType& score_type,
                                                  bool include_unassigned) except + nogil
        void switchScores(ConsensusMap& cmap, Size& counter, bool unassigned_peptides_too) except + nogil
        void switchScores(PeptideIdentificationList& pep_ids, Size& counter) except + nogil
        @staticmethod
        IDScoreSwitcherAlgorithm_IDSwitchResult switchToScoreType(ConsensusMap& cmap,
                                                                  String requested_score_type_as_string,
                                                                  bool include_unassigned) except + nogil
        @staticmethod
        IDScoreSwitcherAlgorithm_IDSwitchResult switchToScoreType(PeptideIdentificationList& pep_ids,
                                                                  String requested_score_type_as_string) except + nogil
        @staticmethod
        void switchBackScoreType(ConsensusMap& cmap,
                                 IDScoreSwitcherAlgorithm_IDSwitchResult isr,
                                 bool include_unassigned) except + nogil
        @staticmethod
        void switchBackScoreType(PeptideIdentificationList& pep_ids,
                                 IDScoreSwitcherAlgorithm_IDSwitchResult isr) except + nogil

cdef extern from "<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>" namespace "OpenMS::IDScoreSwitcherAlgorithm":
    cdef enum IDScoreSwitcherAlgorithm_ScoreType "OpenMS::IDScoreSwitcherAlgorithm::ScoreType":
        # wrap-attach:
        #   IDScoreSwitcherAlgorithm
        RAW
        RAW_EVAL
        PP
        PEP
        FDR
        QVAL

    cdef cppclass IDScoreSwitcherAlgorithm_IDSwitchResult "OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult":
        IDScoreSwitcherAlgorithm_IDSwitchResult() except + nogil
        IDScoreSwitcherAlgorithm_IDSwitchResult(IDScoreSwitcherAlgorithm_IDSwitchResult&) except + nogil
        String original_score_name
        bool original_score_higher_better
        IDScoreSwitcherAlgorithm_ScoreType original_score_type
        bool requested_score_higher_better
        IDScoreSwitcherAlgorithm_ScoreType requested_score_type
        String requested_score_name
        bool score_switched
