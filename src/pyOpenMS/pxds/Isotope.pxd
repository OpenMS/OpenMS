from libcpp cimport bool
from Types cimport *
from String cimport *
from Element cimport *

cdef extern from "<OpenMS/CHEMISTRY/Isotope.h>" namespace "OpenMS":

    cdef cppclass Isotope:

        Isotope() nogil except +
        Isotope(Isotope &) nogil except +

        # detailed constructor
        Isotope(String name,
                String symbol,
                UInt atomic_number,
                UInt neutrons,
                double mono_weight,
                double abundance,
                double half_life,
                DecayMode dm) nogil except +

        # sets unique atomic number
        void setAtomicNumber(UInt atomic_number) nogil except + # wrap-doc:Sets unique atomic number

        # returns the unique atomic number
        UInt getAtomicNumber() nogil except + # wrap-doc:Returns the unique atomic number

        # sets the average weight of the element
        void setAverageWeight(double weight) nogil except + # wrap-doc:Sets the average weight of the element

        # returns the average weight of the element
        double getAverageWeight() nogil except + # wrap-doc:Returns the average weight of the element

        # sets the mono isotopic weight of the element
        void setMonoWeight(double weight) nogil except + # wrap-doc:Sets the mono isotopic weight of the element

        # returns the mono isotopic weight of the element
        double getMonoWeight() nogil except + # wrap-doc:Returns the mono isotopic weight of the element

        # sets the isotope distribution of the element
        void setIsotopeDistribution(IsotopeDistribution isotopes) nogil except + # wrap-doc:Sets the isotope distribution of the element

        # returns the isotope distribution of the element
        IsotopeDistribution getIsotopeDistribution() nogil except + # wrap-doc:Returns the isotope distribution of the element

        # set the name of the element
        void setName(String name) nogil except + # wrap-doc:Sets the name of the element

        # returns the name of the element
        String getName() nogil except + # wrap-doc:Returns the name of the element

        # sets symbol of the element
        void setSymbol(String symbol) nogil except + # wrap-doc:Sets symbol of the element

        # returns symbol of the element
        String getSymbol() nogil except + # wrap-doc:Returns symbol of the element


        const Element* getElement() # wrap-doc:Get corresponding element

        void setHalfLife(double hl) nogil except + # wrap-doc:set isotope half life in seconds
        double getHalfLife() nogil except + # wrap-doc:get isotope half life in seconds
        void setAbundance(double ab) nogil except + # wrap-doc:set isotope natural abundance
        double getAbundance() nogil except + # wrap-doc:get isotope natural abundance
        void setNeutrons(int ne) nogil except + # wrap-doc:set number of neutrons 
        int getNeutrons() nogil except + # wrap-doc:get number of neutrons 
        void setDecayMode(DecayMode dm) nogil except + # wrap-doc:set primary decay mode (for unstable isotopes)
        DecayMode getDecayMode() nogil except + # wrap-doc:get primary decay mode (for unstable isotopes)

        bool isIsotope() nogil except + # wrap-doc:Whether this is an Isotope or an Element (for casting)

        bool isStable() nogil except + # wrap-doc:Whether this is a stable isotope


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
