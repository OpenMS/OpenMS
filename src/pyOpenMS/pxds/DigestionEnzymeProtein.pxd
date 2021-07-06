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
        void setName(String name) nogil except + # wrap-doc:Sets the name of the enzyme

        # returns the name of the Enzyme
        String getName() nogil except + # wrap-doc:Returns the name of the Enzyme

        # sets the synonyms
        void setSynonyms(libcpp_set[String] synonyms) nogil except + # wrap-doc:Sets the synonyms

        # adds a synonym
        void addSynonym(String synonym) nogil except + # wrap-doc:Adds a synonym

        # returns the sysnonyms
        libcpp_set[String] getSynonyms() nogil except + # wrap-doc:Returns the sysnonyms

        # sets the name of the Enzyme as three letter code
        void setRegEx(String three_letter_code) nogil except + # wrap-doc:Sets the name of the Enzyme as three letter code

        # returns the name of the Enzyme as three letter code
        String getRegEx() nogil except + # wrap-doc:Returns the name of the Enzyme as three letter code

        # sets the regex description
        void setRegExDescription(String one_letter_code) nogil except + # wrap-doc:Sets the regex description

        # returns the regex description
        String getRegExDescription() nogil except + # wrap-doc:Returns the regex description

        # sets the N-term gain
        void setNTermGain(EmpiricalFormula value) nogil except + # wrap-doc:Sets the N-term gain

        # sets the C-term gain
        void setCTermGain(EmpiricalFormula value) nogil except + # wrap-doc:Sets the C-term gain

        # returns the N-term gain
        EmpiricalFormula getNTermGain() nogil except + # wrap-doc:Returns the N-term gain

        # returns the C-term gain
        EmpiricalFormula getCTermGain() nogil except + # wrap-doc:Returns the C-term gain

        # sets the PSI ID
        void setPSIID(String value) nogil except + # wrap-doc:Sets the PSI ID

        # returns the PSI ID
        String getPSIID() nogil except + # wrap-doc:Returns the PSI ID

        void setXTandemID(String value) nogil except + # wrap-doc:Sets the X! Tandem enzyme ID
        String getXTandemID() nogil except + # wrap-doc:Returns the X! Tandem enzyme ID

        String getCruxID() nogil except + # wrap-doc:Returns the Crux enzyme ID
        void setCruxID(const String & value) nogil except + # wrap-doc:Sets the Crux enzyme ID

        void setCometID(int value) nogil except + # wrap-doc:Sets the Comet enzyme ID
        int getCometID() nogil except + # wrap-doc:Returns the Comet enzyme ID

        # sets the OMSSA ID
        void setOMSSAID(int value) nogil except + # wrap-doc:Sets the OMSSA ID
        int getOMSSAID() nogil except +

        void setMSGFID(Int value) nogil except +
        Int getMSGFID() nogil except +

        # equality operator
        bool operator==(DigestionEnzymeProtein& Enzyme) nogil except +

        # inequality operator
        bool operator!=(DigestionEnzymeProtein& Enzyme) nogil except +

        # equality operator for cleavage regex
        bool operator==(EmpiricalFormula cleavage_regex) nogil except +

        # equality operator for cleavage regex
        bool operator!=(EmpiricalFormula cleavage_regex) nogil except +

        bool setValueFromFile(const String & key, const String & value) nogil except + # TODO
