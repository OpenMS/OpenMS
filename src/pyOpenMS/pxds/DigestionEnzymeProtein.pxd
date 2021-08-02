from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from Types cimport *
from String cimport *
from EmpiricalFormula cimport *
from DigestionEnzyme cimport *

cdef extern from "<OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>" namespace "OpenMS":

    cdef cppclass DigestionEnzymeProtein(DigestionEnzyme):
        # wrap-inherits:
        #    DigestionEnzyme
        #
        # wrap-doc:
        #   Representation of a digestion enzyme for proteins (protease)

        DigestionEnzymeProtein() nogil except +

        DigestionEnzymeProtein(DigestionEnzymeProtein &) nogil except +

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

        void setCruxID(const String & value) nogil except + # wrap-doc:Sets the Crux enzyme ID
        String getCruxID() nogil except + # wrap-doc:Returns the Crux enzyme ID

        void setCometID(int value) nogil except + # wrap-doc:Sets the Comet enzyme ID
        int getCometID() nogil except + # wrap-doc:Returns the Comet enzyme ID

        void setOMSSAID(int value) nogil except + # wrap-doc:Sets the OMSSA enzyme ID
        int getOMSSAID() nogil except + # wrap-doc:Returns the OMSSA enzyme ID

        void setMSGFID(Int value) nogil except + # wrap-doc:Sets the MSGFPlus enzyme id
        Int getMSGFID() nogil except + # wrap-doc:Returns the MSGFPlus enzyme id

        # equality operator
        bool operator==(DigestionEnzymeProtein& enzyme) nogil except +

        # inequality operator
        bool operator!=(DigestionEnzymeProtein& enzyme) nogil except +

        # order operator
        bool operator<(DigestionEnzymeProtein& enzyme) nogil except +

        # equality operator for cleavage regex
        bool operator==(String cleavage_regex) nogil except +

        # equality operator for cleavage regex
        bool operator!=(String cleavage_regex) nogil except +
