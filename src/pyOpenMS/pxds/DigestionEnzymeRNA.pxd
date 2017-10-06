from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from Types cimport *
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>" namespace "OpenMS":

    cdef cppclass DigestionEnzymeRNA:
        DigestionEnzymeRNA() nogil except + # wrap-ignore

        DigestionEnzymeRNA(DigestionEnzymeRNA) nogil except + # wrap-ignore

        # sets the name of the Enzyme
        void setName(String name) nogil except +

        # returns the name of the Enzyme
        String getName() nogil except +

        # sets the synonyms
        void setSynonyms(libcpp_set[String] synonyms) nogil except +

        # adds a synonym
        void addSynonym(String synonym) nogil except +

        # returns the sysnonyms
        libcpp_set[String] getSynonyms() nogil except +

        # sets the name of the Enzyme as three letter code
        void setRegEx(String three_letter_code) nogil except +

        # returns the name of the Enzyme as three letter code
        String getRegEx() nogil except +

        # sets the regex description
        void setRegExDescription(String one_letter_code) nogil except +

        # returns the regex description
        String getRegExDescription() nogil except +

        # sets the 3' gain
        void setThreePrimeGain(String value) nogil except +

        # sets the 5' gain
        void setFivePrimeGain(String value) nogil except +

        # returns the 3' gain
        String getThreePrimeGain() nogil except +

        # returns the 5' gain
        String getFivePrimeGain() nogil except +

        # equality operator
        bool operator==(DigestionEnzymeRNA& Enzyme) nogil except +

        # inequality operator
        bool operator!=(DigestionEnzymeRNA& Enzyme) nogil except +

        # equality operator for cleavage regex
        bool operator==(EmpiricalFormula cleavage_regex) nogil except +

        # equality operator for cleavage regex
        bool operator!=(EmpiricalFormula cleavage_regex) nogil except +
