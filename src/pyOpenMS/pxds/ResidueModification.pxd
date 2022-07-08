from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS":

    cdef cppclass ResidueModification "OpenMS::ResidueModification":
        # wrap-hash:
        #   getFullId().c_str()

        ResidueModification() nogil except +
        ResidueModification(ResidueModification &) nogil except +

        bool operator==(ResidueModification & modification) nogil except +
        bool operator!=(ResidueModification & modification) nogil except +

        void setId(const String & id_) nogil except + # wrap-doc:Sets the identifier of the modification
        String  getId() nogil except + # wrap-doc:Returns the identifier of the modification
        void setFullId(const String & full_id) nogil except + # wrap-doc:Sets the full identifier (Unimod Accession + origin, if available)
        String  getFullId() nogil except +

        Int getUniModRecordId() nogil except + # wrap-doc:Gets the unimod record id
        void setUniModRecordId(Int id_) nogil except + # wrap-doc:Sets the unimod record id
        String  getUniModAccession() nogil except + # wrap-doc:Returns the unimod accession if available

        void setPSIMODAccession(const String & id_) nogil except + # wrap-doc:Sets the MOD-XXXXX accession of PSI-MOD
        String getPSIMODAccession() nogil except + # wrap-doc:Returns the PSI-MOD accession if available
        void setFullName(const String & full_name) nogil except + # wrap-doc:Sets the full name of the modification; must NOT contain the origin (or . for terminals!)
        String getFullName() nogil except + # wrap-doc:Returns the full name of the modification
        void setName(const String & name) nogil except + # wrap-doc:Sets the name of modification
        String getName() nogil except + # wrap-doc:Returns the PSI-MS-label if available; e.g. Mascot uses this name
        void setTermSpecificity(TermSpecificity term_spec) nogil except + # wrap-doc:Sets the term specificity
        void setTermSpecificity(const String & name) nogil except + # wrap-doc:Sets the terminal specificity using a name

        TermSpecificity getTermSpecificity() nogil except + # wrap-doc:Returns terminal specificity
        String getTermSpecificityName(TermSpecificity ) nogil except + # wrap-doc:Returns the name of the terminal specificity
        void setOrigin(char origin) nogil except + # wrap-doc:Sets the origin (i.e. modified amino acid)
        char  getOrigin() nogil except + # wrap-doc:Returns the origin (i.e. modified amino acid)

        void setSourceClassification(const String & classification) nogil except + # wrap-doc:Classification as defined by the PSI-MOD
        void setSourceClassification(SourceClassification classification) nogil except + # wrap-doc:Sets the source classification
        SourceClassification getSourceClassification() nogil except + # wrap-doc:Returns the source classification, if none was set, it is unspecific
        String getSourceClassificationName(SourceClassification classification) nogil except + # wrap-doc:Returns the classification

        void setAverageMass(double mass) nogil except + # wrap-doc:Sets the average mass
        double getAverageMass() nogil except + # wrap-doc:Returns the average mass if set
        void setMonoMass(double mass) nogil except + # wrap-doc:Sets the monoisotopic mass (this must include the weight of the residue itself!)
        double getMonoMass() nogil except + # wrap-doc:Return the monoisotopic mass, or 0.0 if not set
        void setDiffAverageMass(double mass) nogil except + # wrap-doc:Sets the difference average mass
        double getDiffAverageMass() nogil except + # wrap-doc:Returns the difference average mass, or 0.0 if not set
        void setDiffMonoMass(double mass) nogil except + # wrap-doc:Sets the difference monoisotopic mass
        double getDiffMonoMass() nogil except + # wrap-doc:Returns the diff monoisotopic mass, or 0.0 if not set

        void setFormula(const String & composition) nogil except + # wrap-doc:Sets the formula (no masses will be changed)
        String  getFormula() nogil except + # wrap-doc:Returns the chemical formula if set
        void setDiffFormula(EmpiricalFormula & diff_formula) nogil except + # wrap-doc:Sets diff formula (no masses will be changed)
        EmpiricalFormula  getDiffFormula() nogil except + # wrap-doc:Returns the diff formula if one was set

        void setSynonyms(libcpp_set[ String ] & synonyms) nogil except + # wrap-doc:Sets the synonyms of that modification
        void addSynonym(const String & synonym) nogil except + # wrap-doc:Adds a synonym to the unique list
        libcpp_set[ String ] getSynonyms() nogil except + # wrap-doc:Returns the set of synonyms

        void setNeutralLossDiffFormulas(libcpp_vector[ EmpiricalFormula ] & diff_formulas) nogil except + # wrap-doc:Sets the neutral loss formula
        libcpp_vector[ EmpiricalFormula ] getNeutralLossDiffFormulas() nogil except + # wrap-doc:Returns the neutral loss diff formula (if available)

        void setNeutralLossMonoMasses(libcpp_vector[ double ] mono_masses) nogil except + # wrap-doc:Sets the neutral loss mono weight
        libcpp_vector[ double ] getNeutralLossMonoMasses() nogil except + # wrap-doc:Returns the neutral loss mono weight
        void setNeutralLossAverageMasses(libcpp_vector[ double ] average_masses) nogil except + # wrap-doc:Sets the neutral loss average weight
        libcpp_vector[ double ] getNeutralLossAverageMasses() nogil except + # wrap-doc:Returns the neutral loss average weight

        bool hasNeutralLoss() nogil except + # wrap-doc:Returns true if a neutral loss formula is set
        bool isUserDefined() nogil except + # wrap-doc:Returns true if it is a user-defined modification (empty id)

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS::ResidueModification":
    cdef enum TermSpecificity "OpenMS::ResidueModification::TermSpecificity":
        #wrap-attach:
        #    ResidueModification
        ANYWHERE
        C_TERM
        N_TERM
        PROTEIN_C_TERM
        PROTEIN_N_TERM
        NUMBER_OF_TERM_SPECIFICITY

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS::ResidueModification":
    cdef enum SourceClassification "OpenMS::ResidueModification::SourceClassification":
        #wrap-attach:
        #    ResidueModification
        ARTIFACT
        HYPOTHETICAL
        NATURAL
        POSTTRANSLATIONAL
        MULTIPLE
        CHEMICAL_DERIVATIVE
        ISOTOPIC_LABEL
        PRETRANSLATIONAL
        OTHER_GLYCOSYLATION
        NLINKED_GLYCOSYLATION
        AA_SUBSTITUTION
        OTHER
        NONSTANDARD_RESIDUE
        COTRANSLATIONAL
        OLINKED_GLYCOSYLATION
        UNKNOWN
        NUMBER_OF_SOURCE_CLASSIFICATIONS
