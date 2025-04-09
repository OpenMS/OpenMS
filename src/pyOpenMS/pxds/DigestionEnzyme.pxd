from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from EmpiricalFormula cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzyme.h>" namespace "OpenMS":
    
    cdef cppclass DigestionEnzyme "OpenMS::DigestionEnzyme":
        # wrap-doc:
        #    Base class for digestion enzymes

        DigestionEnzyme(DigestionEnzyme &) except + nogil 

        DigestionEnzyme(const String & name, const String & cleavage_regex, libcpp_set[ String ] & synonyms, String regex_description) except + nogil 

        void setName(const String & name) except + nogil  # wrap-doc:Sets the name of the enzyme

        String getName() except + nogil  # wrap-doc:Returns the name of the enzyme

        void setSynonyms(libcpp_set[ String ] & synonyms) except + nogil  # wrap-doc:Sets the synonyms

        void addSynonym(const String & synonym) except + nogil  # wrap-doc:Adds a synonym

        libcpp_set[ String ] getSynonyms() except + nogil  # wrap-doc:Returns the synonyms

        void setRegEx(const String & cleavage_regex) except + nogil  # wrap-doc:Sets the cleavage regex

        String getRegEx() except + nogil  # wrap-doc:Returns the cleavage regex

        void setRegExDescription(const String & value) except + nogil  # wrap-doc:Sets the regex description

        String getRegExDescription() except + nogil  # wrap-doc:Returns the regex description

        bool operator==(const DigestionEnzyme & enzyme) except + nogil 

        bool operator!=(const DigestionEnzyme & enzyme) except + nogil 

        bool operator<(const DigestionEnzyme & enzyme) except + nogil 

        bool setValueFromFile(String key, String value) except + nogil  # wrap-doc:Sets the value of a member variable based on an entry from an input file
