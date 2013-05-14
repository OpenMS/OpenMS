from libcpp.string cimport string as libcpp_string
from libcpp.set cimport set as libcpp_set
#from InstrumentSettings cimport *
from CVTermList cimport *
from Peak1D cimport *
from Map cimport *

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS":

    cdef cppclass Precursor(Peak1D):

        # wrap-inherits:
        #    Peak1D

        #####################################################################
        # we do not declare inheritance from CVTermList, because cython
        # has problems inheriting overloaded methods (in this case
        # replaceCVTerms).
        # Instead we declare the derived methods manually:

        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                             String accession
                            ) nogil except +


        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                            ) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms()
        void addCVTerm(CVTerm & term)                   nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        #####################################################################
        # we do not declare inheritance from MetaInfoInterface, because cython
        # has problems inheriting overloaded methods (in this case
        # e.g. getKeys and getMetaValue).
        # Instead we declare the derived methods manually:

        void getKeys(libcpp_vector[String] & keys)
        void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +


        #####################################################################

        Precursor()           nogil except +
        Precursor(Precursor)           nogil except +

        # returns a mutable reference to the activation methods
        libcpp_set[ActivationMethod] getActivationMethods() nogil except +
        # sets the activation methods
        void setActivationMethods(libcpp_set[ActivationMethod] activation_methods) nogil except +

        double getActivationEnergy() nogil except +
        void setActivationEnergy(double activation_energy) nogil except +

        double getIsolationWindowLowerOffset() nogil except +
        void setIsolationWindowLowerOffset(double bound) nogil except +

        double getIsolationWindowUpperOffset() nogil except +
        void setIsolationWindowUpperOffset(double bound) nogil except +

        int getCharge() nogil except +
        void setCharge(int charge) nogil except +

        #Non-mutable access to possible charge states
        libcpp_vector[int] getPossibleChargeStates() nogil except +
        #Sets the possible charge states
        void setPossibleChargeStates(libcpp_vector[int] possible_charge_states) nogil except +

        # Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged
        DoubleReal getUnchargedMass() nogil except +

        bool operator==(Precursor)  nogil except +
        bool operator!=(Precursor)  nogil except +


cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS::Precursor":
    cdef enum ActivationMethod:
      CID,                      #< Collision-induced dissociation
      PSD,                      #< Post-source decay
      PD,                       #< Plasma desorption
      SID,                      #< Surface-induced dissociation
      BIRD,                             #< Blackbody infrared radiative dissociation
      ECD,                              #< Electron capture dissociation
      IMD,                              #< Infrared multiphoton dissociation
      SORI,                             #< Sustained off-resonance irradiation
      HCID,                             #< High-energy collision-induced dissociation
      LCID,                             #< Low-energy collision-induced dissociation
      PHD,                              #< Photodissociation
      ETD,                              #< Electron transfer dissociation
      PQD,                              #< Pulsed q dissociation
      SIZE_OF_ACTIVATIONMETHOD

