from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from Types cimport *
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>" namespace "OpenMS":

    cdef cppclass DigestionEnzymeProtein:
        DigestionEnzymeProtein() nogil except + # wrap-ignore

        DigestionEnzymeProtein(DigestionEnzymeProtein) nogil except + # wrap-ignore

        # detailed constructor
        DigestionEnzymeProtein(String name,
                               String cleavage_regex,
                               libcpp_set[String] synonyms,
                               String regex_description,
                               EmpiricalFormula n_term_gain,
                               EmpiricalFormula c_term_gain,
                               String psi_id,
                               String xtandem_id,
                               UInt comet_id,
                               UInt omssa_id) nogil except +

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

        # sets the N-term gain
        void setNTermGain(EmpiricalFormula value) nogil except +

        # sets the C-term gain
        void setCTermGain(EmpiricalFormula value) nogil except +

        # returns the N-term gain
        EmpiricalFormula getNTermGain() nogil except +

        # returns the C-term gain
        EmpiricalFormula getCTermGain() nogil except +

        # sets the PSI ID
        void setPSIID(String value) nogil except +

        # returns the PSI ID
        String getPSIID() nogil except +

        void setXTandemID(String value) nogil except +
        String getXTandemID() nogil except +

        # sets the Comet ID
        void setCometID(int value) nogil except +

        # returns the Comet ID
        int getCometID() nogil except +

        # sets the OMSSA ID
        void setOMSSAID(int value) nogil except +

        # returns the OMSSA ID
        int getOMSSAID() nogil except +

        # equality operator
        bool operator==(DigestionEnzymeProtein& Enzyme) nogil except +

        # inequality operator
        bool operator!=(DigestionEnzymeProtein& Enzyme) nogil except +

        # equality operator for cleavage regex
        bool operator==(EmpiricalFormula cleavage_regex) nogil except +

        # equality operator for cleavage regex
        bool operator!=(EmpiricalFormula cleavage_regex) nogil except +
