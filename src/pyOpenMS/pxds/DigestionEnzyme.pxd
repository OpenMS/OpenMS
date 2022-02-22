from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from EmpiricalFormula cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzyme.h>" namespace "OpenMS":
    
    cdef cppclass DigestionEnzyme "OpenMS::DigestionEnzyme":
        # wrap-doc:
        #     Base class for digestion enzymes

        DigestionEnzyme(DigestionEnzyme &) nogil except +

        DigestionEnzyme(const String & name, const String & cleavage_regex, libcpp_set[ String ] & synonyms, String regex_description) nogil except +

        void setName(const String & name) nogil except + # wrap-doc:Sets the name of the enzyme

        String getName() nogil except + # wrap-doc:Returns the name of the enzyme

        void setSynonyms(libcpp_set[ String ] & synonyms) nogil except + # wrap-doc:Sets the synonyms

        void addSynonym(const String & synonym) nogil except + # wrap-doc:Adds a synonym

        libcpp_set[ String ] getSynonyms() nogil except + # wrap-doc:Returns the synonyms

        void setRegEx(const String & cleavage_regex) nogil except + # wrap-doc:Sets the cleavage regex

        String getRegEx() nogil except + # wrap-doc:Returns the cleavage regex

        void setRegExDescription(const String & value) nogil except + # wrap-doc:Sets the regex description

        String getRegExDescription() nogil except + # wrap-doc:Returns the regex description

        bool operator==(const DigestionEnzyme & enzyme) nogil except +

        bool operator!=(const DigestionEnzyme & enzyme) nogil except +

        bool operator<(const DigestionEnzyme & enzyme) nogil except +

        bool setValueFromFile(String key, String value) nogil except + # wrap-doc:Sets the value of a member variable based on an entry from an input file
