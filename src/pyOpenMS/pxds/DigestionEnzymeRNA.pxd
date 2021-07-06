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

        bool setValueFromFile(const String & key, const String & value) nogil except +

        # sets the name of the Enzyme
        void setName(String name) nogil except + # wrap-doc:Sets the name of the enzyme

        # returns the name of the Enzyme
        String getName() nogil except + # wrap-doc:Returns the name of the enzyme

        # sets the synonyms
        void setSynonyms(libcpp_set[String] synonyms) nogil except + # wrap-doc:Sets the synonyms

        # adds a synonym
        void addSynonym(String synonym) nogil except + # wrap-doc:Adds a synonym

        # returns the sysnonyms
        libcpp_set[String] getSynonyms() nogil except + # wrap-doc:Returns the synonyms

        # sets the name of the Enzyme as three letter code 
        void setRegEx(String three_letter_code) nogil except + # wrap-doc:Sets the name of the enzyme as three letter code

        # returns the name of the Enzyme as three letter code
        String getRegEx() nogil except + # wrap-doc:Returns the name of the enzyme as three letter code

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
        bool operator==(DigestionEnzymeRNA& Enzyme) nogil except + 

        # inequality operator
        bool operator!=(DigestionEnzymeRNA& Enzyme) nogil except + 

        # equality operator for cleavage regex
        bool operator==(EmpiricalFormula cleavage_regex) nogil except + # wrap-doc:Equality operator for cleavage regex

        # equality operator for cleavage regex
        bool operator!=(EmpiricalFormula cleavage_regex) nogil except + # wrap-doc:Equality operator for cleavage regex
