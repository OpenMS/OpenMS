from libcpp.map cimport map as libcpp_map
from Types cimport *
from CVTermList cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS":

    cdef cppclass Precursor(Peak1D, CVTermList):
        # wrap-inherits:
        #   Peak1D
        #   CVTermList
        # wrap-doc:
        #  Precursor meta information
        #  
        #  This class contains precursor information:
        #
        #  - isolation window
        #  - activation
        #  - selected ion (m/z, intensity, charge, possible charge states)
        #  - ion mobility drift time

        Precursor() except + nogil 
        Precursor(Precursor &) except + nogil 

        libcpp_set[ActivationMethod] getActivationMethods() except + nogil  # wrap-doc:Returns the activation methods
        void setActivationMethods(libcpp_set[ActivationMethod] activation_methods) except + nogil  # wrap-doc:Sets the activation methods

        double getActivationEnergy() except + nogil  # wrap-doc:Returns the activation energy (in electronvolt)
        void setActivationEnergy(double activation_energy) except + nogil  # wrap-doc:Sets the activation energy (in electronvolt)

        double getIsolationWindowLowerOffset() except + nogil  # wrap-doc:Returns the lower offset from the target m/z
        void setIsolationWindowLowerOffset(double bound) except + nogil  # wrap-doc:Sets the lower offset from the target m/z

        double getDriftTime() except + nogil  # wrap-doc:Returns the ion mobility drift time in milliseconds (-1 means it is not set)
        void setDriftTime(double drift_time) except + nogil  # wrap-doc:Sets the ion mobility drift time in milliseconds

        double getIsolationWindowUpperOffset() except + nogil  # wrap-doc:Returns the upper offset from the target m/z
        void setIsolationWindowUpperOffset(double bound) except + nogil  # wrap-doc:Sets the upper offset from the target m/z

        double getDriftTimeWindowLowerOffset() except + nogil  # wrap-doc:Returns the lower offset from the target ion mobility in milliseconds
        void setDriftTimeWindowLowerOffset(double drift_time) except + nogil  # wrap-doc:Sets the lower offset from the target ion mobility
        double getDriftTimeWindowUpperOffset() except + nogil  # wrap-doc:Returns the upper offset from the target ion mobility in milliseconds
        void setDriftTimeWindowUpperOffset(double drift_time) except + nogil  # wrap-doc:Sets the upper offset from the target ion mobility

        int getCharge() except + nogil  # wrap-doc:Returns the charge
        void setCharge(int charge) except + nogil  # wrap-doc:Sets the charge

        libcpp_vector[int] getPossibleChargeStates() except + nogil  # wrap-doc:Returns the possible charge states
        void setPossibleChargeStates(libcpp_vector[int] possible_charge_states) except + nogil  # wrap-doc:Sets the possible charge states

        double getUnchargedMass() except + nogil  # wrap-doc:Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged

        bool operator==(Precursor)  except + nogil 
        bool operator!=(Precursor)  except + nogil 

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS::Precursor":
    cdef enum ActivationMethod:
      # wrap-doc:
      #  Enum for activation/fragmentation methods of mass spectra
      #  
      #  - CID: Collision-induced dissociation (MS:1000133) (also CAD; parent term, but unless otherwise stated often used as synonym for trap-type CID)
      #  - PSD: Post-source decay
      #  - PD: Plasma desorption
      #  - SID: Surface-induced dissociation
      #  - BIRD: Blackbody infrared radiative dissociation
      #  - ECD: Electron capture dissociation (MS:1000250)
      #  - IMD: Infrared multiphoton dissociation
      #  - SORI: Sustained off-resonance irradiation
      #  - HCID: High-energy collision-induced dissociation
      #  - LCID: Low-energy collision-induced dissociation
      #  - PHD: Photodissociation
      #  - ETD: Electron transfer dissociation
      #  - ETciD: Electron transfer and collision-induced dissociation (MS:1003182)
      #  - EThcD: Electron transfer and higher-energy collision dissociation (MS:1002631)
      #  - PQD: Pulsed q dissociation (MS:1000599)
      #  - TRAP: trap-type collision-induced dissociation (MS:1002472)
      #  - HCD: beam-type collision-induced dissociation (MS:1000422)
      #  - INSOURCE: in-source collision-induced dissociation (MS:1001880)
      #  - LIFT: Bruker proprietary method (MS:1002000)
      CID,
      PSD,
      PD,
      SID,
      BIRD,
      ECD,
      IMD,
      SORI,
      HCID,
      LCID,
      PHD,
      ETD,
      ETciD,
      EThcD,
      PQD,
      TRAP,
      HCD,
      INSOURCE,
      LIFT,
      SIZE_OF_ACTIVATIONMETHOD
