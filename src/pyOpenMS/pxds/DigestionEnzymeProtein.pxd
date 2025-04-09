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
        #   DigestionEnzyme
        #
        # wrap-doc:
        #  Representation of a digestion enzyme for proteins (protease)

        DigestionEnzymeProtein() except + nogil 

        DigestionEnzymeProtein(DigestionEnzymeProtein &) except + nogil 

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
                               UInt omssa_id) except + nogil 

        # sets the N-term gain
        void setNTermGain(EmpiricalFormula value) except + nogil  # wrap-doc:Sets the N-term gain

        # sets the C-term gain
        void setCTermGain(EmpiricalFormula value) except + nogil  # wrap-doc:Sets the C-term gain

        # returns the N-term gain
        EmpiricalFormula getNTermGain() except + nogil  # wrap-doc:Returns the N-term gain

        # returns the C-term gain
        EmpiricalFormula getCTermGain() except + nogil  # wrap-doc:Returns the C-term gain

        # sets the PSI ID
        void setPSIID(String value) except + nogil  # wrap-doc:Sets the PSI ID

        # returns the PSI ID
        String getPSIID() except + nogil  # wrap-doc:Returns the PSI ID

        void setXTandemID(String value) except + nogil  # wrap-doc:Sets the X! Tandem enzyme ID
        String getXTandemID() except + nogil  # wrap-doc:Returns the X! Tandem enzyme ID

        void setCometID(int value) except + nogil  # wrap-doc:Sets the Comet enzyme ID
        int getCometID() except + nogil  # wrap-doc:Returns the Comet enzyme ID

        void setOMSSAID(int value) except + nogil  # wrap-doc:Sets the OMSSA enzyme ID
        int getOMSSAID() except + nogil  # wrap-doc:Returns the OMSSA enzyme ID

        void setMSGFID(Int value) except + nogil  # wrap-doc:Sets the MSGFPlus enzyme id
        Int getMSGFID() except + nogil  # wrap-doc:Returns the MSGFPlus enzyme id

        # equality operator
        bool operator==(DigestionEnzymeProtein& enzyme) except + nogil 

        # inequality operator
        bool operator!=(DigestionEnzymeProtein& enzyme) except + nogil 

        # order operator
        bool operator<(DigestionEnzymeProtein& enzyme) except + nogil 

        # equality operator for cleavage regex
        bool operator==(String cleavage_regex) except + nogil 

        # equality operator for cleavage regex
        bool operator!=(String cleavage_regex) except + nogil 
