from libcpp cimport bool
from Types cimport *
from String cimport *
from Element cimport *

cdef extern from "<OpenMS/CHEMISTRY/Isotope.h>" namespace "OpenMS":

    cdef cppclass Isotope:

        Isotope() except + nogil
        Isotope(Isotope &) except + nogil

        # detailed constructor
        Isotope(String name,
                String symbol,
                UInt atomic_number,
                UInt neutrons,
                double mono_weight,
                double abundance,
                double half_life,
                DecayMode dm) except + nogil

        # sets unique atomic number
        void setAtomicNumber(UInt atomic_number) except + nogil # wrap-doc:Sets unique atomic number

        # returns the unique atomic number
        UInt getAtomicNumber() except + nogil # wrap-doc:Returns the unique atomic number

        # sets the average weight of the element
        void setAverageWeight(double weight) except + nogil # wrap-doc:Sets the average weight of the element

        # returns the average weight of the element
        double getAverageWeight() except + nogil # wrap-doc:Returns the average weight of the element

        # sets the mono isotopic weight of the element
        void setMonoWeight(double weight) except + nogil # wrap-doc:Sets the mono isotopic weight of the element

        # returns the mono isotopic weight of the element
        double getMonoWeight() except + nogil # wrap-doc:Returns the mono isotopic weight of the element

        # sets the isotope distribution of the element
        void setIsotopeDistribution(IsotopeDistribution isotopes) except + nogil # wrap-doc:Sets the isotope distribution of the element

        # returns the isotope distribution of the element
        IsotopeDistribution getIsotopeDistribution() except + nogil # wrap-doc:Returns the isotope distribution of the element

        # set the name of the element
        void setName(String name) except + nogil # wrap-doc:Sets the name of the element

        # returns the name of the element
        String getName() except + nogil # wrap-doc:Returns the name of the element

        # sets symbol of the element
        void setSymbol(String symbol) except + nogil # wrap-doc:Sets symbol of the element

        # returns symbol of the element
        String getSymbol() except + nogil # wrap-doc:Returns symbol of the element


        const Element* getElement() except + nogil # wrap-doc:Get corresponding element

        void setHalfLife(double hl) except + nogil # wrap-doc:set isotope half life in seconds
        double getHalfLife() except + nogil # wrap-doc:get isotope half life in seconds
        void setAbundance(double ab) except + nogil # wrap-doc:set isotope natural abundance
        double getAbundance() except + nogil # wrap-doc:get isotope natural abundance
        void setNeutrons(int ne) except + nogil # wrap-doc:set number of neutrons 
        int getNeutrons() except + nogil # wrap-doc:get number of neutrons 
        void setDecayMode(DecayMode dm) except + nogil # wrap-doc:set primary decay mode (for unstable isotopes)
        DecayMode getDecayMode() except + nogil # wrap-doc:get primary decay mode (for unstable isotopes)

        bool isIsotope() except + nogil # wrap-doc:Whether this is an Isotope or an Element (for casting)

        bool isStable() except + nogil # wrap-doc:Whether this is a stable isotope


cdef extern from "<OpenMS/CHEMISTRY/Isotope.h>" namespace "OpenMS::Isotope":

    cdef enum DecayMode:
      # wrap-attach:
      #   Isotope
      NONE = 0,     # No decay (stable isotope)
      UNKNOWN,      # Unknown / Unspecified decay mode
      ALPHA,        # Alpha decay
      BETA_PLUS,    # Beta plus decay
      BETA_MINUS,   # Beta minus decay
      PROTON,       # Proton emission
      SIZE_OF_DECAYMODE
