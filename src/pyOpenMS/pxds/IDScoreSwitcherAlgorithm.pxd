from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from DefaultParamHandler cimport *
from PeptideIdentification cimport *
from PeptideIdentificationList cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>" namespace "OpenMS":

    cdef cppclass IDScoreSwitcherAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        #
        # wrap-doc:
        #  Switches identification scores between different score types.
        #
        #  This class provides functionality to switch the main scoring type used in peptide
        #  or protein identification data. It supports switching between different score types,
        #  such as raw scores, E-values, posterior probabilities, posterior error probabilities,
        #  FDR, and q-values.
        #
        #  The class can detect score types in both main scores and metavalues, making it useful
        #  for finding PEP (Posterior Error Probability) values regardless of where they are stored.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    from pyopenms import *
        #    idsa = IDScoreSwitcherAlgorithm()
        #
        #    # Check if a score name is of PEP type
        #    score_type = pep_id.getScoreType()
        #    is_pep = idsa.isScoreType(score_type, ScoreType.PEP)
        #
        #    # Check if higher score is better for a score type
        #    higher_better = idsa.isScoreTypeHigherBetter(ScoreType.PEP)  # Returns False
        #
        #    # Get all known score names
        #    names = idsa.getScoreNames()

        IDScoreSwitcherAlgorithm() except + nogil
        IDScoreSwitcherAlgorithm(IDScoreSwitcherAlgorithm &) except + nogil  # wrap-ignore

        bool isScoreType(const String & score_name, ScoreType type_) except + nogil
            # wrap-doc:
            #  Checks if the given score name corresponds to a specific score type.
            #
            #  Performs a case-insensitive comparison and handles the "_score" suffix.
            #
            #  :param score_name: The name of the score to check
            #  :param type_: The ScoreType to compare against
            #  :return: True if the score name matches the given ScoreType

        bool isScoreTypeHigherBetter(ScoreType score_type) except + nogil
            # wrap-doc:
            #  Determines whether a higher score is better for the given ScoreType.
            #
            #  :param score_type: The score type to check
            #  :return: True if higher scores are better, False otherwise

        libcpp_vector[String] getScoreNames() except + nogil
            # wrap-doc:
            #  Gets a vector of all score names that are used in OpenMS.
            #
            #  :return: Vector of score names (e.g., "q-value", "ln(hyperscore)")



cdef extern from "<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>" namespace "OpenMS::IDScoreSwitcherAlgorithm":

    cdef enum ScoreType "OpenMS::IDScoreSwitcherAlgorithm::ScoreType":
        # wrap-doc:
        #  Score type hierarchy for MS identification scores.
        #
        #  Used to categorize and detect different types of scores:
        #  - RAW: Raw search engine scores (e.g., hyperscore)
        #  - RAW_EVAL: Raw scores with E-value (e.g., expect score)
        #  - PP: Posterior Probability
        #  - PEP: Posterior Error Probability
        #  - FDR: False Discovery Rate
        #  - QVAL: Q-value

        # wrap-attach:
        #  IDScoreSwitcherAlgorithm

        RAW  # wrap-doc:Raw score (e.g., hyperscore, XTandem)
        RAW_EVAL  # wrap-doc:Raw score with E-value (e.g., expect score)
        PP  # wrap-doc:Posterior Probability
        PEP  # wrap-doc:Posterior Error Probability
        FDR  # wrap-doc:False Discovery Rate
        QVAL  # wrap-doc:Q-value

    ScoreType toScoreTypeEnum(String score_type) except + nogil
        # wrap-attach:
        #  IDScoreSwitcherAlgorithm
        # wrap-doc:
        #  Converts a string representation of a score type to ScoreType enum.
        #
        #  :param score_type: String like "PEP", "q-value", "FDR", etc.
        #  :return: The corresponding ScoreType enum value
        #  :raises: Exception if the string doesn't match any known score type
