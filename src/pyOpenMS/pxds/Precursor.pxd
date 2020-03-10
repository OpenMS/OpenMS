from Types cimport *
from CVTermList cimport *
from Peak1D cimport *
from Map cimport *

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS":

    cdef cppclass Precursor(Peak1D, CVTermList):
        # wrap-inherits:
        #    Peak1D
        #    CVTermList

        Precursor() nogil except +
        Precursor(Precursor) nogil except +

        # returns a mutable reference to the activation methods
        libcpp_set[ActivationMethod] getActivationMethods() nogil except +
        # sets the activation methods
        void setActivationMethods(libcpp_set[ActivationMethod] activation_methods) nogil except +

        double getActivationEnergy() nogil except +
        void setActivationEnergy(double activation_energy) nogil except +

        double getIsolationWindowLowerOffset() nogil except +
        void setIsolationWindowLowerOffset(double bound) nogil except +

        double getDriftTime() nogil except +
        void setDriftTime(double drift_time) nogil except +

        double getIsolationWindowUpperOffset() nogil except +
        void setIsolationWindowUpperOffset(double bound) nogil except +

        double getDriftTimeWindowLowerOffset() nogil except +
        void setDriftTimeWindowLowerOffset(double drift_time) nogil except +
        double getDriftTimeWindowUpperOffset() nogil except +
        void setDriftTimeWindowUpperOffset(double drift_time) nogil except +

        int getCharge() nogil except +
        void setCharge(int charge) nogil except +

        #Non-mutable access to possible charge states
        libcpp_vector[int] getPossibleChargeStates() nogil except +
        #Sets the possible charge states
        void setPossibleChargeStates(libcpp_vector[int] possible_charge_states) nogil except +

        # Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged
        double getUnchargedMass() nogil except +

        bool operator==(Precursor)  nogil except +
        bool operator!=(Precursor)  nogil except +

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS::Precursor":
    cdef enum ActivationMethod:
      CID,                      #< Collision-induced dissociation
      PSD,                      #< Post-source decay
      PD,                       #< Plasma desorption
      SID,                      #< Surface-induced dissociation
      BIRD,                     #< Blackbody infrared radiative dissociation
      ECD,                      #< Electron capture dissociation
      IMD,                      #< Infrared multiphoton dissociation
      SORI,                     #< Sustained off-resonance irradiation
      HCID,                     #< High-energy collision-induced dissociation
      LCID,                     #< Low-energy collision-induced dissociation
      PHD,                      #< Photodissociation
      ETD,                      #< Electron transfer dissociation
      PQD,                      #< Pulsed q dissociation
      SIZE_OF_ACTIVATIONMETHOD

