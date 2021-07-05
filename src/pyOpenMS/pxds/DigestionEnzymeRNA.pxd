from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from Types cimport *
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>" namespace "OpenMS":

    cdef cppclass DigestionEnzymeRNA:
        # wrap-inherits:
        #    DigestionEnzyme
        #
        # wrap-doc:
        #   Representation of a digestion enzyme for RNA (RNase)
        #   -----
        #   The cutting sites of these enzymes are defined using two different mechanisms:
        #   First, a single regular expression that is applied to strings of unmodified RNA sequence and defines cutting sites via zero-length matches (using lookahead/lookbehind assertions).
        #   This is the same mechanism that is used for proteases (@see @ref ProteaseDigestion).
        #   However, due to the complex notation involved, this approach is not practical for modification-aware digestion.
        #   Thus, the second mechanism uses two regular expressions ("cuts after"/"cuts before"), which are applied to the short codes (e.g. "m6A") of sequential ribonucleotides.
        #   If both expressions match, then there is a cutting site between the two ribonucleotides.
        #   -----
        #   There is support for terminal (5'/3') modifications that may be generated on fragments as a result of RNase cleavage.
        #   A typical example is 3'-phosphate, resulting from cleavage of the phosphate backbone.

        DigestionEnzymeRNA() nogil except +

        DigestionEnzymeRNA(DigestionEnzymeRNA) nogil except +

        bool setValueFromFile(const String & key, const String & value) nogil except + # wrap-doc:Set the value of a member variable based on an entry from an input file

        # sets the name of the Enzyme
        void setName(String name) nogil except + # wrap-doc:Sets the name of the Enzyme

        # returns the name of the Enzyme
        String getName() nogil except + # wrap-doc:Returns the name of the Enzyme

        # sets the synonyms
        void setSynonyms(libcpp_set[String] synonyms) nogil except + # wrap-doc:Sets the synonyms

        # adds a synonym
        void addSynonym(String synonym) nogil except + # wrap-doc:Adds a synonym

        # returns the synonyms
        libcpp_set[String] getSynonyms() nogil except + # wrap-doc:Returns the synonyms

        # sets the name of the Enzyme as three letter code
        void setRegEx(String three_letter_code) nogil except + # wrap-doc:Sets the name of the Enzyme as three letter code

        # returns the name of the Enzyme as three letter code
        String getRegEx() nogil except + # wrap-doc:Returns the name of the Enzyme as three letter code

        # sets the regex description
        void setRegExDescription(String one_letter_code) nogil except + # wrap-doc:Sets the regex description

        # returns the regex description
        String getRegExDescription() nogil except + # wrap-doc:Returns the regex description

        # sets the 3' gain
        void setThreePrimeGain(String value) nogil except + # wrap-doc:Sets the 3' gain

        # sets the 5' gain
        void setFivePrimeGain(String value) nogil except + # wrap-doc:Sets the 5' gain

        # returns the 3' gain
        String getThreePrimeGain() nogil except + # wrap-doc:Returns the 3' gain

        # returns the 5' gain
        String getFivePrimeGain() nogil except + # wrap-doc:Returns the 5' gain

        # equality operator
        bool operator==(DigestionEnzymeRNA& Enzyme) nogil except + # wrap-doc:Equality operator

        # inequality operator
        bool operator!=(DigestionEnzymeRNA& Enzyme) nogil except + # wrap-doc:Inequality operator

        # equality operator for cleavage regex
        bool operator==(EmpiricalFormula cleavage_regex) nogil except + # wrap-doc:Equality operator for cleavage regex

        # equality operator for cleavage regex
        bool operator!=(EmpiricalFormula cleavage_regex) nogil except + # wrap-doc:Equality operator for cleavage regex

