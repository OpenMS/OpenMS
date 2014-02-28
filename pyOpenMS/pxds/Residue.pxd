from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from Types cimport *
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/Residue.h>" namespace "OpenMS":

    cdef cppclass Residue:

        Residue() nogil except +
        Residue(Residue) nogil except + # wrap-ignore

        # detailed constructor
        Residue(String name,
                String three_letter_code,
                String one_letter_code,
                EmpiricalFormula formula) nogil except +

        # Internal
        EmpiricalFormula getInternalToFull() nogil except +

        DoubleReal getInternalToFullAverageWeight() nogil except +

        DoubleReal getInternalToFullMonoWeight() nogil except +

        # N-terminal
        EmpiricalFormula getNTerminalToFull() nogil except +

        DoubleReal getNTerminalToFullAverageWeight() nogil except +

        DoubleReal getNTerminalToFullMonoWeight() nogil except +

        # C-terminal
        EmpiricalFormula getCTerminalToFull() nogil except +

        DoubleReal getCTerminalToFullAverageWeight() nogil except +

        DoubleReal getCTerminalToFullMonoWeight() nogil except +

        # b ion
        EmpiricalFormula getBIonToFull() nogil except +

        DoubleReal getBIonToFullAverageWeight() nogil except +

        DoubleReal getBIonToFullMonoWeight() nogil except +

        # a ion
        EmpiricalFormula getAIonToFull() nogil except +

        DoubleReal getAIonToFullAverageWeight() nogil except +

        DoubleReal getAIonToFullMonoWeight() nogil except +

        # y ion
        EmpiricalFormula getYIonToFull() nogil except +

        DoubleReal getYIonToFullAverageWeight() nogil except +

        DoubleReal getYIonToFullMonoWeight() nogil except +

        # c ion
        EmpiricalFormula getCIonToFull() nogil except +

        DoubleReal getCIonToFullAverageWeight() nogil except +

        DoubleReal getCIonToFullMonoWeight() nogil except +

        # c-1 ion
        EmpiricalFormula getCIonMinusOneToFull() nogil except +

        DoubleReal getCIonMinusOneToFullAverageWeight() nogil except +

        DoubleReal getCIonMinusOneToFullMonoWeight() nogil except +

        # c+1 ion
        EmpiricalFormula getCIonPlusOneToFull() nogil except +

        DoubleReal getCIonPlusOneToFullAverageWeight() nogil except +

        DoubleReal getCIonPlusOneToFullMonoWeight() nogil except +

        # c+2 ion
        EmpiricalFormula getCIonPlusTwoToFull() nogil except +

        DoubleReal getCIonPlusTwoToFullAverageWeight() nogil except +

        DoubleReal getCIonPlusTwoToFullMonoWeight() nogil except +

        # x ion
        EmpiricalFormula getXIonToFull() nogil except +

        DoubleReal getXIonToFullAverageWeight() nogil except +

        DoubleReal getXIonToFullMonoWeight() nogil except +

        # z ion
        EmpiricalFormula getZIonToFull() nogil except +

        DoubleReal getZIonToFullAverageWeight() nogil except +

        DoubleReal getZIonToFullMonoWeight() nogil except +

        # z-1 ion
        EmpiricalFormula getZIonMinusOneToFull() nogil except +

        DoubleReal getZIonMinusOneToFullAverageWeight() nogil except +

        DoubleReal getZIonMinusOneToFullMonoWeight() nogil except +

        # z+1 ion
        EmpiricalFormula getZIonPlusOneToFull() nogil except +

        DoubleReal getZIonPlusOneToFullAverageWeight() nogil except +

        DoubleReal getZIonPlusOneToFullMonoWeight() nogil except +

        # z+2 ion
        EmpiricalFormula getZIonPlusTwoToFull() nogil except +

        DoubleReal getZIonPlusTwoToFullAverageWeight() nogil except +

        DoubleReal getZIonPlusTwoToFullMonoWeight() nogil except +

        # returns the ion name given as a residue type
        String getResidueTypeName(ResidueType res_type) nogil except +

        # sets the name of the residue
        void setName(String name) nogil except +

        # returns the name of the residue
        String getName() nogil except +

        # sets the short name of the residue, this name is used in the PeptideSequence for output
        void setShortName(String short_name) nogil except +

        # returns the short name of the residue
        String getShortName() nogil except +

        # sets the synonyms
        void setSynonyms(libcpp_set[String] synonyms) nogil except +

        # adds a synonym
        void addSynonym(String synonym) nogil except +

        # returns the sysnonyms
        libcpp_set[String] getSynonyms() nogil except +

        # sets the name of the residue as three letter code
        void setThreeLetterCode(String three_letter_code) nogil except +

        # returns the name of the residue as three letter code
        String getThreeLetterCode() nogil except +

        # sets the name as one letter code
        void setOneLetterCode(String one_letter_code) nogil except +

        # returns the name as one letter code
        String getOneLetterCode() nogil except +

        # adds a neutral loss formula
        void addLossFormula(EmpiricalFormula) nogil except +

        # sets the neutral loss formulas
        void setLossFormulas(libcpp_vector[EmpiricalFormula]) nogil except +

        # adds N-terminal losses
        void addNTermLossFormula(EmpiricalFormula) nogil except +

        # sets the N-terminal losses
        void setNTermLossFormulas(libcpp_vector[EmpiricalFormula]) nogil except +

        # returns the neutral loss formulas
        libcpp_vector[EmpiricalFormula] getLossFormulas() nogil except +

        # returns N-terminal loss formulas
        libcpp_vector[EmpiricalFormula] getNTermLossFormulas() nogil except +

        # set the neutral loss molecule name
        void setLossNames(libcpp_vector[String] name) nogil except +

        # sets the N-terminal loss names
        void setNTermLossNames(libcpp_vector[String] name) nogil except +

        # add neutral loss molecule name
        void addLossName(String name) nogil except +

        # adds a N-terminal loss name
        void addNTermLossName(String name) nogil except +

        # gets neutral loss name (if there is one, else returns an empty string)
        libcpp_vector[String] getLossNames() nogil except +

        # returns the N-terminal loss names
        libcpp_vector[String] getNTermLossNames() nogil except +

        # set empirical formula of the residue (must be full, with N and C-terminus)
        void setFormula(EmpiricalFormula formula) nogil except +

        # returns the empirical formula of the residue
        EmpiricalFormula getFormula(ResidueType res_type) nogil except +

        # sets average weight of the residue (must be full, with N and C-terminus)
        void setAverageWeight(DoubleReal weight) nogil except +

        # returns average weight of the residue
        DoubleReal getAverageWeight(ResidueType res_type) nogil except +

        # sets mono weight of the residue (must be full, with N and C-terminus)
        void setMonoWeight(DoubleReal weight) nogil except +

        # returns mono weight of the residue
        DoubleReal getMonoWeight(ResidueType res_type) nogil except +

        # sets by the name, this mod should be present in ModificationsDB
        void setModification(String name) nogil except +

        # returns the name of the modification to the modification
        String getModification() nogil except +

        # sets the low mass marker ions as a vector of formulas
        void setLowMassIons(libcpp_vector[EmpiricalFormula] low_mass_ions) nogil except +

        # returns a vector of formulas with the low mass markers of the residue
        libcpp_vector[EmpiricalFormula] getLowMassIons() nogil except +

        # sets the residue sets the amino acid is contained in
        void setResidueSets(libcpp_set[String] residues_sets) nogil except +

        # adds a residue set to the residue sets
        void addResidueSet(String residue_sets) nogil except +

        # returns the residue sets this residue is contained in
        libcpp_set[String] getResidueSets() nogil except +

        # true if the residue has neutral loss
        bool hasNeutralLoss() nogil except +

        # true if N-terminal neutral losses are set
        bool hasNTermNeutralLosses() nogil except +

        # equality operator
        bool operator==(Residue & residue) nogil except +

        # inequality operator
        bool operator!=(Residue & residue) nogil except +

        # equality operator for one letter code
        bool operator==(char one_letter_code) nogil except +

        # equality operator for one letter code
        bool operator!=(char one_letter_code) nogil except +

        # returns the pka of the residue
        DoubleReal getPka() nogil except +

        # returns the pkb of the residue
        DoubleReal getPkb() nogil except +

        # returns the pkc of the residue if it exists otherwise -1
        DoubleReal getPkc() nogil except +

        # calculates the isoelectric point using the pk* values
        DoubleReal getPiValue() nogil except +

        # sets the pka of the residue
        void setPka(DoubleReal value) nogil except +

        # sets the pkb of the residue
        void setPkb(DoubleReal value) nogil except +

        # sets the pkc of the residue
        void setPkc(DoubleReal value) nogil except +

        # returns the side chain basicity
        DoubleReal getSideChainBasicity() nogil except +

        # sets the side chain basicity
        void setSideChainBasicity(DoubleReal gb_sc) nogil except +

        # returns the backbone basicitiy if located in N-terminal direction
        DoubleReal getBackboneBasicityLeft() nogil except +

        # sets the N-terminal direction backbone basicitiy
        void setBackboneBasicityLeft(DoubleReal gb_bb_l) nogil except +

        # returns the C-terminal direction backbone basicitiy
        DoubleReal getBackboneBasicityRight() nogil except +

        # sets the C-terminal direction backbone basicity
        void setBackboneBasicityRight(DoubleReal gb_bb_r) nogil except +

        # true if the residue is a modified one
        bool isModified() nogil except +

        # true if the residue is contained in the set
        bool isInResidueSet(String residue_set) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/Residue.h>" namespace "OpenMS::Residue":

    cdef enum ResidueType:
      # wrap-attach:
      #   Residue
      Full = 0,           # with N-terminus and C-terminus
      Internal,           # internal, without any termini
      NTerminal,           # only N-terminus
      CTerminal,           # only C-terminus
      AIon,           # N-terminus up to the C-alpha/carbonyl carbon bond
      BIon,           # N-terminus up to the peptide bond
      CIonMinusOne,           # N-terminus up to the amide/C-alpha bond
      CIon,           # N-terminus up to the amide/C-alpha bond
      CIonPlusOne,           # N-terminus up to the amide/C-alpha bond
      CIonPlusTwo,           # N-terminus up to the amide/C-alpha bond
      XIon,           # amide/C-alpha bond up to the C-terminus
      YIon,           # peptide bond up to the C-terminus
      ZIonMinusOne,           # C-alpha/carbonyl carbon bond
      ZIon,            # C-alpha/carbonyl carbon bond
      ZIonPlusOne,            # C-alpha/carbonyl carbon bond
      ZIonPlusTwo,            # C-alpha/carbonyl carbon bond
      SizeOfResidueType

