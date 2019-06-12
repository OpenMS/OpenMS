from Types cimport *
from String cimport *
from EmpiricalFormula cimport *


cdef extern from "<OpenMS/CHEMISTRY/Ribonucleotide.h>" namespace "OpenMS::Ribonucleotide":

    cdef enum TermSpecificityNuc:
      # wrap-attach:
      # Ribonucleotide
      ANYWHERE = 0,
      FIVE_PRIME,
      THREE_PRIME,
      NUMBER_OF_TERM_SPECIFICITY


cdef extern from "<OpenMS/CHEMISTRY/Ribonucleotide.h>" namespace "OpenMS":

    cdef cppclass Ribonucleotide:
        # wrap-hash:
        #   getName().c_str()

        Ribonucleotide() nogil except +
        Ribonucleotide(Ribonucleotide) nogil except + # wrap-ignore

        # detailed constructor
        Ribonucleotide(String name,
                String code,
                String new_code,
                String html_code,
                EmpiricalFormula formula,
                char origin,
                double mono_mass,
                double avg_mass,
                TermSpecificityNuc term_spec,
                EmpiricalFormula baseloss_formula) nogil except +

 
        # return the short name
        String getCode() nogil except +

        # set the short name
        void setCode(String code) nogil except +

        # sets the name of the ribonucleotide
        void setName(String name) nogil except +

        # returns the name of the ribonucleotide
        String getName() nogil except +

        # set empirical formula of the ribonucleotide (must be full, with N and C-terminus)
        void setFormula(EmpiricalFormula formula) nogil except +

        # returns the empirical formula of the residue
        EmpiricalFormula getFormula() nogil except +

        # sets average mass of the ribonucleotide
        void setAvgMass(double avg_mass) nogil except +

        # returns average mass of the ribonucleotide
        double getAvgMass() nogil except +

        # sets monoisotopic mass of the ribonucleotide
        void setMonoMass(double mono_mass) nogil except +

        # returns monoisotopic mass of the ribonucleotide
        double getMonoMass() nogil except +

        # return the new code
        String getNewCode() nogil except +

        # set the new code
        void setNewCode(String code) nogil except +

        # Get the code of the unmodified base (e.g., "A", "C", ...)
        char getOrigin() nogil except +

        # Set the code of the unmodified base (e.g., "A", "C", ...)
        void setOrigin(char origin) nogil except +

        # Set the HTML (RNAMods) code
        void setHTMLCode(String html_code) nogil except +

        # Get the HTML (RNAMods) code
        String getHTMLCode() nogil except +

        # Set the terminal specificity
        void setTermSpecificity(TermSpecificityNuc term_spec) nogil except +

        # Get the terminal specificity
        TermSpecificityNuc getTermSpecificity() nogil except +

        # Get sum formula after loss of the nucleobase
        EmpiricalFormula getBaselossFormula() nogil except +

        # Set sum formula after loss of the nucleobase
        void setBaselossFormula(EmpiricalFormula formula) nogil except +

        # equality operator
        bool operator==(Ribonucleotide & ribonucleotide) nogil except +

        # true if the ribonucleotide is a modified one
        bool isModified() nogil except +
