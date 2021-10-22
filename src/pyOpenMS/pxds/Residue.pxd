from Types cimport *
from String cimport *
from EmpiricalFormula cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/Residue.h>" namespace "OpenMS":

    cdef cppclass Residue:
        # wrap-hash:
        #   getName().c_str()

        Residue() nogil except +
        Residue(Residue &) nogil except +

        # detailed constructor
        Residue(String name,
                String three_letter_code,
                String one_letter_code,
                EmpiricalFormula formula) nogil except +

        # Conversions
        EmpiricalFormula getInternalToFull() nogil except +
        EmpiricalFormula getInternalToNTerm() nogil except +
        EmpiricalFormula getInternalToCTerm() nogil except +
        EmpiricalFormula getInternalToAIon() nogil except +
        EmpiricalFormula getInternalToBIon() nogil except +
        EmpiricalFormula getInternalToCIon() nogil except +
        EmpiricalFormula getInternalToXIon() nogil except +
        EmpiricalFormula getInternalToYIon() nogil except +
        EmpiricalFormula getInternalToZIon() nogil except +

        String getResidueTypeName(ResidueType res_type) nogil except + # wrap-doc:Returns the ion name given as a residue type

        void setName(String name) nogil except + # wrap-doc:Sets the name of the residue

        String getName() nogil except + # wrap-doc:Returns the name of the residue

        void setSynonyms(libcpp_set[String] synonyms) nogil except + # wrap-doc:Sets the synonyms

        void addSynonym(String synonym) nogil except + # wrap-doc:Adds a synonym

        libcpp_set[String] getSynonyms() nogil except + # wrap-doc:Returns the sysnonyms

        void setThreeLetterCode(String three_letter_code) nogil except + # wrap-doc:Sets the name of the residue as three letter code

        String getThreeLetterCode() nogil except + # wrap-doc:Returns the name of the residue as three letter code

        void setOneLetterCode(String one_letter_code) nogil except + # wrap-doc:Sets the name as one letter code

        String getOneLetterCode() nogil except + # wrap-doc:Returns the name as one letter code

        void addLossFormula(EmpiricalFormula) nogil except + # wrap-doc:Adds a neutral loss formula

        void setLossFormulas(libcpp_vector[EmpiricalFormula]) nogil except + # wrap-doc:Sets the neutral loss formulas

        void addNTermLossFormula(EmpiricalFormula) nogil except + # wrap-doc:Adds N-terminal losses

        void setNTermLossFormulas(libcpp_vector[EmpiricalFormula]) nogil except + # wrap-doc:Sets the N-terminal losses

        libcpp_vector[EmpiricalFormula] getLossFormulas() nogil except + # wrap-doc:Returns the neutral loss formulas

        libcpp_vector[EmpiricalFormula] getNTermLossFormulas() nogil except + # wrap-doc:Returns N-terminal loss formulas

        void setLossNames(libcpp_vector[String] name) nogil except + # wrap-doc:Sets the neutral loss molecule name

        void setNTermLossNames(libcpp_vector[String] name) nogil except + # wrap-doc:Sets the N-terminal loss names

        void addLossName(String name) nogil except + # wrap-doc:Adds neutral loss molecule name

        void addNTermLossName(String name) nogil except + # wrap-doc:Adds a N-terminal loss name

        libcpp_vector[String] getLossNames() nogil except + # wrap-doc:Gets neutral loss name (if there is one, else returns an empty string)

        libcpp_vector[String] getNTermLossNames() nogil except + # wrap-doc:Returns the N-terminal loss names

        void setFormula(EmpiricalFormula formula) nogil except + # wrap-doc:Sets empirical formula of the residue (must be full, with N and C-terminus)

        EmpiricalFormula getFormula() nogil except + # wrap-doc:Returns the empirical formula of the residue
        EmpiricalFormula getFormula(ResidueType res_type) nogil except +

        void setAverageWeight(double weight) nogil except + # wrap-doc:Sets average weight of the residue (must be full, with N and C-terminus)

        double getAverageWeight() nogil except + # wrap-doc:Returns average weight of the residue
        double getAverageWeight(ResidueType res_type) nogil except +

        void setMonoWeight(double weight) nogil except + # wrap-doc:Sets monoisotopic weight of the residue (must be full, with N and C-terminus)

        double getMonoWeight() nogil except + # wrap-doc:Returns monoisotopic weight of the residue
        double getMonoWeight(ResidueType res_type) nogil except +

        const ResidueModification * getModification() nogil except +

        # setModification by pointer is not here since a copy would be made whose memory is not handled by the ModificationDB
        void setModification(String name) nogil except + # wrap-doc:Sets the modification by name; the mod should be present in ModificationsDB

        void setModification(const ResidueModification& mod) nogil except + # wrap-doc:Sets the modification by a ResidueModification object; checks if present in ModificationsDB and adds if not.

        void setModificationByDiffMonoMass(double diffMonoMass) nogil except + # wrap-doc:Sets the modification by monoisotopic mass difference in Da; checks if present in ModificationsDB with tolerance and adds a "user-defined" modification if not (for later lookups).

        String getModificationName() nogil except + # wrap-doc:Returns the name of the modification to the modification

        void setLowMassIons(libcpp_vector[EmpiricalFormula] low_mass_ions) nogil except + # wrap-doc:Sets the low mass marker ions as a vector of formulas

        libcpp_vector[EmpiricalFormula] getLowMassIons() nogil except + # wrap-doc:Returns a vector of formulas with the low mass markers of the residue

        void setResidueSets(libcpp_set[String] residues_sets) nogil except + # wrap-doc:Sets the residue sets the amino acid is contained in

        void addResidueSet(String residue_sets) nogil except + # wrap-doc:Adds a residue set to the residue sets

        libcpp_set[String] getResidueSets() nogil except + # wrap-doc:Returns the residue sets this residue is contained in

        bool hasNeutralLoss() nogil except + # wrap-doc:True if the residue has neutral loss

        bool hasNTermNeutralLosses() nogil except + # wrap-doc:True if N-terminal neutral losses are set

        # equality operator
        bool operator==(Residue & residue) nogil except +

        # inequality operator
        bool operator!=(Residue & residue) nogil except +

        # equality operator for one letter code
        bool operator==(char one_letter_code) nogil except +

        # equality operator for one letter code
        bool operator!=(char one_letter_code) nogil except +

        double getPka() nogil except + # wrap-doc:Returns the pka of the residue

        double getPkb() nogil except + # wrap-doc:Returns the pkb of the residue

        double getPkc() nogil except + # wrap-doc:Returns the pkc of the residue if it exists otherwise -1

        double getPiValue() nogil except + # wrap-doc:Calculates the isoelectric point using the pk values

        void setPka(double value) nogil except + # wrap-doc:Sets the pka of the residue

        void setPkb(double value) nogil except + # wrap-doc:Sets the pkb of the residue

        void setPkc(double value) nogil except + # wrap-doc:Sets the pkc of the residue

        double getSideChainBasicity() nogil except + # wrap-doc:Returns the side chain basicity

        void setSideChainBasicity(double gb_sc) nogil except + # wrap-doc:Sets the side chain basicity

        double getBackboneBasicityLeft() nogil except + # wrap-doc:Returns the backbone basicitiy if located in N-terminal direction

        void setBackboneBasicityLeft(double gb_bb_l) nogil except + # wrap-doc:Sets the N-terminal direction backbone basicitiy

        double getBackboneBasicityRight() nogil except + # wrap-doc:Returns the C-terminal direction backbone basicitiy

        void setBackboneBasicityRight(double gb_bb_r) nogil except + # wrap-doc:Sets the C-terminal direction backbone basicity

        bool isModified() nogil except + # wrap-doc:True if the residue is a modified one

        bool isInResidueSet(String residue_set) nogil except + # wrap-doc:True if the residue is contained in the set

        char residueTypeToIonLetter(ResidueType res_type) nogil except + # wrap-doc:Helper for mapping residue types to letters for Text annotations and labels

cdef extern from "<OpenMS/CHEMISTRY/Residue.h>" namespace "OpenMS::Residue":

    cdef enum ResidueType:
      # wrap-attach:
      #   Residue
      Full = 0,       # with N-terminus and C-terminus
      Internal,       # internal, without any termini
      NTerminal,      # only N-terminus
      CTerminal,      # only C-terminus
      AIon,           # MS:1001229 N-terminus up to the C-alpha/carbonyl carbon bond
      BIon,           # MS:1001224 N-terminus up to the peptide bond
      CIon,           # MS:1001231 N-terminus up to the amide/C-alpha bond
      XIon,           # MS:1001228 amide/C-alpha bond up to the C-terminus
      YIon,           # MS:1001220 peptide bond up to the C-terminus
      ZIon,           # MS:1001230 C-alpha/carbonyl carbon bond
      Precursor_ion,  # MS:1001523 Precursor ion
      BIonMinusH20,   # MS:1001222 b ion without water
      YIonMinusH20,   # MS:1001223 y ion without water
      BIonMinusNH3,   # MS:1001232 b ion without ammonia
      YIonMinusNH3,   # MS:1001233 y ion without ammonia
      NonIdentified,  # MS:1001240 Non-identified ion
      Unannotated,    # no stored annotation
      SizeOfResidueType
