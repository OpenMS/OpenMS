from Types cimport *
from CVTermList cimport *
from Peak1D cimport *
from Map cimport *

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS":

    cdef cppclass Precursor(Peak1D, CVTermList):
        # wrap-inherits:
        #    Peak1D
        #    CVTermList
        # wrap-doc:
        #   Precursor meta information
        #   -----
        #   This class contains precursor information:
        #
        #   - isolation window
        #   - activation
        #   - selected ion (m/z, intensity, charge, possible charge states)
        #   - ion mobility drift time

        Precursor() nogil except +
        Precursor(Precursor &) nogil except +

        libcpp_set[ActivationMethod] getActivationMethods() nogil except + # wrap-doc:Returns the activation methods
        void setActivationMethods(libcpp_set[ActivationMethod] activation_methods) nogil except + # wrap-doc:Sets the activation methods

        double getActivationEnergy() nogil except + # wrap-doc:Returns the activation energy (in electronvolt)
        void setActivationEnergy(double activation_energy) nogil except + # wrap-doc:Sets the activation energy (in electronvolt)

        double getIsolationWindowLowerOffset() nogil except + # wrap-doc:Returns the lower offset from the target m/z
        void setIsolationWindowLowerOffset(double bound) nogil except + # wrap-doc:Sets the lower offset from the target m/z

        double getDriftTime() nogil except + # wrap-doc:Returns the ion mobility drift time in milliseconds (-1 means it is not set)
        void setDriftTime(double drift_time) nogil except + # wrap-doc:Sets the ion mobility drift time in milliseconds

        double getIsolationWindowUpperOffset() nogil except + # wrap-doc:Returns the upper offset from the target m/z
        void setIsolationWindowUpperOffset(double bound) nogil except + # wrap-doc:Sets the upper offset from the target m/z

        double getDriftTimeWindowLowerOffset() nogil except + # wrap-doc:Returns the lower offset from the target ion mobility in milliseconds
        void setDriftTimeWindowLowerOffset(double drift_time) nogil except + # wrap-doc:Sets the lower offset from the target ion mobility
        double getDriftTimeWindowUpperOffset() nogil except + # wrap-doc:Returns the upper offset from the target ion mobility in milliseconds
        void setDriftTimeWindowUpperOffset(double drift_time) nogil except + # wrap-doc:Sets the upper offset from the target ion mobility

        int getCharge() nogil except + # wrap-doc:Returns the charge
        void setCharge(int charge) nogil except + # wrap-doc:Sets the charge

        libcpp_vector[int] getPossibleChargeStates() nogil except + # wrap-doc:Returns the possible charge states
        void setPossibleChargeStates(libcpp_vector[int] possible_charge_states) nogil except + # wrap-doc:Sets the possible charge states

        double getUnchargedMass() nogil except + # wrap-doc:Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged

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
      TRAP,                     #< trap-type collision-induced dissociation (MS:1002472)
      HCD,                      #< beam-type collision-induced dissociation (MS:1000422) / HCD
      INSOURCE,                 #< in-source collision-induced dissociation (MS:1001880)
      LIFT,                     #< Bruker proprietary method (MS:1002000)
      SIZE_OF_ACTIVATIONMETHOD

