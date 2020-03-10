from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from EmpiricalFormula cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzyme.h>" namespace "OpenMS":
    
    cdef cppclass DigestionEnzyme "OpenMS::DigestionEnzyme":
        DigestionEnzyme(DigestionEnzyme) nogil except +
        DigestionEnzyme(const String & name, const String & cleavage_regex, libcpp_set[ String ] & synonyms, String regex_description) nogil except +
        void setName(const String & name) nogil except +
        String getName() nogil except +
        void setSynonyms(libcpp_set[ String ] & synonyms) nogil except +
        void addSynonym(const String & synonym) nogil except +
        libcpp_set[ String ] getSynonyms() nogil except +
        void setRegEx(const String & cleavage_regex) nogil except +
        String getRegEx() nogil except +
        void setRegExDescription(const String & value) nogil except +
        String getRegExDescription() nogil except +
        bool operator==(const DigestionEnzyme & enzyme) nogil except +
        bool operator!=(const DigestionEnzyme & enzyme) nogil except +
        bool operator==(const String & cleavage_regex) nogil except +
        bool operator!=(const String & cleavage_regex) nogil except +
        bool operator<(const DigestionEnzyme & enzyme) nogil except +
        bool setValueFromFile(const String & key, const String & value) nogil except +

