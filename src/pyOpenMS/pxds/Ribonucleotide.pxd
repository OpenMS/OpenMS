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
        #  getName().c_str()

        Ribonucleotide() except + nogil 
        Ribonucleotide(Ribonucleotide &) except + nogil 

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
                EmpiricalFormula baseloss_formula) except + nogil 

 
        String getCode() except + nogil  # wrap-doc:Returns the short name

        void setCode(String code) except + nogil  # wrap-doc:Sets the short name

        void setName(String name) except + nogil  # wrap-doc:Sets the name of the ribonucleotide

        String getName() except + nogil  # wrap-doc:Returns the name of the ribonucleotide

        void setFormula(EmpiricalFormula formula) except + nogil  # wrap-doc:Sets empirical formula of the ribonucleotide (must be full, with N and C-terminus)

        EmpiricalFormula getFormula() except + nogil  # wrap-doc:Returns the empirical formula of the residue

        void setAvgMass(double avg_mass) except + nogil  # wrap-doc:Sets average mass of the ribonucleotide

        double getAvgMass() except + nogil  # wrap-doc:Returns average mass of the ribonucleotide

        void setMonoMass(double mono_mass) except + nogil  # wrap-doc:Sets monoisotopic mass of the ribonucleotide

        double getMonoMass() except + nogil  # wrap-doc:Returns monoisotopic mass of the ribonucleotide

        String getNewCode() except + nogil  # wrap-doc:Returns the new code

        void setNewCode(String code) except + nogil  # wrap-doc:Sets the new code

        char getOrigin() except + nogil  # wrap-doc:Returns the code of the unmodified base (e.g., "A", "C", ...)

        void setOrigin(char origin) except + nogil  # wrap-doc:Sets the code of the unmodified base (e.g., "A", "C", ...)

        void setHTMLCode(String html_code) except + nogil  # wrap-doc:Sets the HTML (RNAMods) code

        String getHTMLCode() except + nogil  # wrap-doc:Returns the HTML (RNAMods) code

        void setTermSpecificity(TermSpecificityNuc term_spec) except + nogil  # wrap-doc:Sets the terminal specificity

        TermSpecificityNuc getTermSpecificity() except + nogil  # wrap-doc:Returns the terminal specificity

        EmpiricalFormula getBaselossFormula() except + nogil  # wrap-doc:Returns sum formula after loss of the nucleobase

        void setBaselossFormula(EmpiricalFormula formula) except + nogil  # wrap-doc:Sets sum formula after loss of the nucleobase

        bool operator==(Ribonucleotide & ribonucleotide) except + nogil  # wrap-doc:Equality operator

        bool isModified() except + nogil  # wrap-doc:True if the ribonucleotide is a modified one
