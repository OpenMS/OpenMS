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
                double average_weight,
                double abundance,
                IsotopeDistribution isotopes) nogil except +

        const Element* getElement() # wrap-doc: Get corresponding element

        void setHalfLife(double hl) nogil except + # wrap-doc: set isotope half life in seconds
        double getHalfLife() nogil except + # wrap-doc: get isotope half life in seconds
        void setAbundance(double ab) nogil except + # wrap-doc: set isotope natural abundance
        double getAbundance() nogil except + # wrap-doc: get isotope natural abundance
        void setNeutrons(int ne) nogil except + # wrap-doc: set number of neutrons 
        int getNeutrons() nogil except + # wrap-doc: get number of neutrons 
        void setDecayMode(DecayMode dm) nogil except + # wrap-doc: set primary decay mode (for unstable isotopes)
        DecayMode getDecayMode() nogil except + # wrap-doc: get primary decay mode (for unstable isotopes)

        virtual bool isIsotope() nogil except + # wrap-doc: Whether this is an Isotope or an Element (for casting)

        bool isStable() nogil except + # wrap-doc: Whether this is a stable isotope


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
