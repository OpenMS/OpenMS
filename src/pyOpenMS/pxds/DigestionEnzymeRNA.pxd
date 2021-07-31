from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from Types cimport *
from String cimport *
from DigestionEnzyme cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>" namespace "OpenMS":

    cdef cppclass DigestionEnzymeRNA(DigestionEnzyme):
        # wrap-inherits:
        #    DigestionEnzyme
        #
        # wrap-doc:
        #   Representation of a digestion enzyme for RNA (RNase)
        #   -----
        #   The cutting sites of these enzymes are defined using two different mechanisms:
        #   First, a single regular expression that is applied to strings of unmodified RNA sequence and defines cutting sites via zero-length matches (using lookahead/lookbehind assertions).
        #   This is the same mechanism that is used for proteases (see ProteaseDigestion).
        #   However, due to the complex notation involved, this approach is not practical for modification-aware digestion.
        #   Thus, the second mechanism uses two regular expressions ("cuts after"/"cuts before"), which are applied to the short codes (e.g. "m6A") of sequential ribonucleotides.
        #   If both expressions match, then there is a cutting site between the two ribonucleotides.
        #   -----
        #   There is support for terminal (5'/3') modifications that may be generated on fragments as a result of RNase cleavage.
        #   A typical example is 3'-phosphate, resulting from cleavage of the phosphate backbone.

        DigestionEnzymeRNA() nogil except +

        DigestionEnzymeRNA(DigestionEnzymeRNA &) nogil except +

        # sets the "cuts after ..." regular expression
        void setCutsAfterRegEx(String value) nogil except + # wrap-doc:Sets the "cuts after ..." regular expression

        # returns the "cuts after ..." regular expression
        String getCutsAfterRegEx() nogil except + # wrap-doc:Returns the "cuts after ..." regular expression

        # sets the "cuts before ..." regular expression
        void setCutsBeforeRegEx(String value) nogil except + # wrap-doc:Sets the "cuts before ..." regular expression

        # returns the "cuts before ..." regular expression
        String getCutsBeforeRegEx() nogil except + # wrap-doc:Returns the "cuts before ..." regular expression

        # sets the 3' gain
        void setThreePrimeGain(String value) nogil except + # wrap-doc:Sets the 3' gain

        # sets the 5' gain
        void setFivePrimeGain(String value) nogil except + # wrap-doc:Sets the 5' gain

        # returns the 3' gain
        String getThreePrimeGain() nogil except + # wrap-doc:Returns the 3' gain

        # returns the 5' gain
        String getFivePrimeGain() nogil except + # wrap-doc:Returns the 5' gain
