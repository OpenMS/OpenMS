from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from String cimport *


# Define the enum first so it can be used by the Scores class
cdef extern from "<OpenMS/ANALYSIS/ID/Scores.h>" namespace "OpenMS::Scores":

    cdef enum class IDScoreType "OpenMS::Scores::Type":
        # wrap-attach:
        #    Scores
        # wrap-doc:
        #  Hierarchy of possible score types in MS identification
        #
        #  - RAW: Raw score, e.g., search engine specific scores like hyperscore
        #  - RAW_EVAL: Raw score with E-value, e.g., expect score
        #  - PP: Posterior probability
        #  - PEP: Posterior error probability
        #  - FDR: False discovery rate
        #  - QVAL: Q-value
        RAW,
        RAW_EVAL,
        PP,
        PEP,
        FDR,
        QVAL


cdef extern from "<OpenMS/ANALYSIS/ID/Scores.h>" namespace "OpenMS":

    cdef cppclass Scores:
        # wrap-doc:
        #  Utility class for score type handling in identification workflows.
        #
        #  This class provides centralized handling of score types used in peptide and protein
        #  identification. It defines the hierarchy of score types (raw scores, E-values,
        #  posterior probabilities, etc.) and provides utility methods for score type
        #  conversion, comparison, and lookup.

        Scores() except + nogil
        Scores(Scores &) except + nogil

        @staticmethod
        bool isScoreType(const String& score_name, IDScoreType type) except + nogil
            # wrap-doc:
            #  Checks if the given score name corresponds to a specific score type
            #
            #  :param score_name: The name of the score to check
            #  :param type: The IDScoreType to compare against
            #  :returns: True if the score name matches the given IDScoreType

        @staticmethod
        IDScoreType toType(const String& score_type) except + nogil
            # wrap-doc:
            #  Converts a string representation of a score type to an IDScoreType enum
            #
            #  :param score_type: The string representation of the score type
            #  :returns: The corresponding IDScoreType enum value
            #  :raises: Exception::MissingInformation if the score_type string is not recognized

        @staticmethod
        bool isHigherBetter(IDScoreType type) except + nogil
            # wrap-doc:
            #  Determines whether a higher score is better for the given score type
            #
            #  :param type: The score type to check
            #  :returns: True if a higher score is better

        @staticmethod
        libcpp_vector[String] getScoreNames() except + nogil
            # wrap-doc:
            #  Gets a vector of all score names that are used in OpenMS
            #
            #  :returns: A vector of all score names (e.g., "q-value", "ln(hyperscore)")

        # Note: getNamesForType and findTypeByName are not wrapped as they use std::set
        # which is more complex to wrap and less useful in Python
