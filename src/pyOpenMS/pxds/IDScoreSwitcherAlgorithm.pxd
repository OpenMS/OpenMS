from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from DefaultParamHandler cimport *
from PeptideIdentification cimport *
from PeptideIdentificationList cimport *
from ConsensusMap cimport *
from Scores cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>" namespace "OpenMS":

    cdef cppclass IDScoreSwitcherAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #  Algorithm to switch identification scores within identification or consensus feature maps
        #
        #  This class provides functionality to switch the main scoring type used in peptide or protein
        #  identification data. It supports switching between different score types, such as raw scores,
        #  E-values, posterior probabilities, posterior error probabilities, FDR, and q-values.

        IDScoreSwitcherAlgorithm() except + nogil
        IDScoreSwitcherAlgorithm(IDScoreSwitcherAlgorithm &) except + nogil

        bool isScoreType(const String& score_name, IDScoreType type) except + nogil
            # wrap-doc:
            #  Checks if the given score name corresponds to a specific score type
            #
            #  :param score_name: The name of the score to check
            #  :param type: The IDScoreType to compare against
            #  :returns: True if the score name matches the given IDScoreType

        ScoreSearchResult findScoreType(PeptideIdentification& id, IDScoreType score_type) except + nogil
            # wrap-doc:
            #  Searches for a score type in a PeptideIdentification
            #
            #  Returns a ScoreSearchResult indicating whether the main score is of the
            #  requested type, and if not, searches for scores of that type in the
            #  meta values of the first hit.
            #
            #  :param id: The PeptideIdentification to analyze
            #  :param score_type: The IDScoreType to search for (e.g., IDScoreType.PEP)
            #  :returns: ScoreSearchResult with is_main_score_type and score_name fields

        @staticmethod
        IDScoreType toScoreTypeEnum(String score_type) except + nogil
            # wrap-doc:
            #  Converts a string representation of a score type to an IDScoreType enum
            #
            #  :param score_type: The string representation of the score type
            #  :returns: The corresponding IDScoreType enum value
            #  :raises: Exception::MissingInformation if the score_type string is not recognized

        bool isScoreTypeHigherBetter(IDScoreType score_type) except + nogil
            # wrap-doc:
            #  Determines whether a higher score type is better given an IDScoreType enum
            #
            #  :param score_type: The score type to check
            #  :returns: True if a higher score type is better

        libcpp_vector[String] getScoreNames() except + nogil
            # wrap-doc:
            #  Gets a vector of all score names that are used in OpenMS
            #
            #  :returns: A vector of all score names (e.g., "q-value", "ln(hyperscore)")

        void switchToGeneralScoreType(PeptideIdentificationList& pep_ids, IDScoreType type, Size& counter) except + nogil
            # wrap-doc:
            #  Switches the score type of a PeptideIdentificationList to a general score type
            #
            #  :param pep_ids: The PeptideIdentificationList whose scores need to be switched
            #  :param type: The desired general score type to switch to
            #  :param counter: A reference to a counter that will be incremented for each peptide identification processed

        void switchToGeneralScoreType(ConsensusMap& cmap, IDScoreType type, Size& counter, bool unassigned_peptides_too) except + nogil
            # wrap-doc:
            #  Switches the score type of a ConsensusMap to a general score type
            #
            #  :param cmap: The ConsensusMap containing peptide identifications whose scores need to be switched
            #  :param type: The desired general score type to switch to
            #  :param counter: A reference to a counter that will be incremented for each peptide identification processed
            #  :param unassigned_peptides_too: Whether to include unassigned peptides in the score switching process

        void switchScores(PeptideIdentificationList& pep_ids, Size& counter) except + nogil
            # wrap-doc:
            #  Switches the scores of peptide identifications
            #
            #  :param pep_ids: The peptide identifications whose scores need to be switched
            #  :param counter: A reference to a counter that will be incremented for each peptide identification processed

        void switchScores(ConsensusMap& cmap, Size& counter, bool unassigned_peptides_too) except + nogil
            # wrap-doc:
            #  Switches the scores of peptide identifications in a ConsensusMap
            #
            #  :param cmap: The ConsensusMap containing peptide identifications whose scores need to be switched
            #  :param counter: A reference to a counter that will be incremented for each peptide identification processed
            #  :param unassigned_peptides_too: Whether to include unassigned peptides in the score switching process

        @staticmethod
        IDSwitchResult switchToScoreType(PeptideIdentificationList& pep_ids, String requested_score_type_as_string) except + nogil
            # wrap-doc:
            #  Switches the score type of peptide identifications to the requested type
            #
            #  :param pep_ids: A vector of PeptideIdentification objects to be processed
            #  :param requested_score_type_as_string: The desired score type as a string (e.g., "RAW", "PEP", "q-value")
            #  :returns: IDSwitchResult containing details about the original and requested score types

        @staticmethod
        IDSwitchResult switchToScoreType(ConsensusMap& cmap, String requested_score_type_as_string, bool include_unassigned) except + nogil
            # wrap-doc:
            #  Switches the score type of a ConsensusMap to the requested score type
            #
            #  :param cmap: The ConsensusMap object whose score types are to be switched
            #  :param requested_score_type_as_string: The desired score type as a string
            #  :param include_unassigned: Whether to include unassigned IDs in the score switch
            #  :returns: IDSwitchResult containing information about the score switch operation

        @staticmethod
        void switchBackScoreType(PeptideIdentificationList& pep_ids, IDSwitchResult isr) except + nogil
            # wrap-doc:
            #  Reverts the scoring type of peptide identifications to their original scores
            #
            #  :param pep_ids: A vector of PeptideIdentification objects to be updated
            #  :param isr: An IDSwitchResult object containing information about the score switch state

        @staticmethod
        void switchBackScoreType(ConsensusMap& cmap, IDSwitchResult isr, bool include_unassigned) except + nogil
            # wrap-doc:
            #  Reverts the score type of a ConsensusMap to its original type
            #
            #  :param cmap: The ConsensusMap object whose scores will be modified
            #  :param isr: The IDSwitchResult containing information about the score switch
            #  :param include_unassigned: Whether to include unassigned PSMs in the score switching process


cdef extern from "<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>" namespace "OpenMS::IDScoreSwitcherAlgorithm":

    cdef cppclass ScoreSearchResult "OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult":
        # wrap-doc:
        #  Structure to hold score detection results for any IDScoreType
        #
        #  Used by findScoreType() to return whether the main score is of a
        #  requested type and the name of the score (either main score name
        #  or metavalue name if found there).

        ScoreSearchResult() except + nogil
        ScoreSearchResult(ScoreSearchResult &) except + nogil

        bool is_main_score_type  # True if main score is already of the requested type
        String score_name        # Name of score to use (main score name or metavalue name)

    cdef cppclass IDSwitchResult "OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult":
        # wrap-doc:
        #  Structure holding score switching information
        #
        #  Contains both the original and requested score details, including
        #  score names, their orientation (whether higher scores are better),
        #  and score types before and after the switch.

        IDSwitchResult() except + nogil
        IDSwitchResult(IDSwitchResult &) except + nogil

        String original_score_name
        bool original_score_higher_better
        IDScoreType original_score_type
        bool requested_score_higher_better
        IDScoreType requested_score_type
        String requested_score_name
        bool score_switched
