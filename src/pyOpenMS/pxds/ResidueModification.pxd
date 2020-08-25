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
        ResidueModification(ResidueModification) nogil except +

        bool operator==(ResidueModification & modification) nogil except +
        bool operator!=(ResidueModification & modification) nogil except +

        void setId(const String & id_) nogil except +
        String  getId() nogil except +
        void setFullId(const String & full_id) nogil except +
        String  getFullId() nogil except +

        Int getUniModRecordId() nogil except +
        void setUniModRecordId(Int id_) nogil except +
        String  getUniModAccession() nogil except +

        void setPSIMODAccession(const String & id_) nogil except +
        String getPSIMODAccession() nogil except +
        void setFullName(const String & full_name) nogil except +
        String getFullName() nogil except +
        void setName(const String & name) nogil except +
        String getName() nogil except +
        void setTermSpecificity(TermSpecificity term_spec) nogil except +
        void setTermSpecificity(const String & name) nogil except +

        TermSpecificity getTermSpecificity() nogil except +
        String getTermSpecificityName(TermSpecificity ) nogil except +
        void setOrigin(char origin) nogil except +
        char  getOrigin() nogil except +

        void setSourceClassification(const String & classification) nogil except +
        void setSourceClassification(SourceClassification classification) nogil except +
        SourceClassification getSourceClassification() nogil except +
        String getSourceClassificationName(SourceClassification classification) nogil except +

        void setAverageMass(double mass) nogil except +
        double getAverageMass() nogil except +
        void setMonoMass(double mass) nogil except +
        double getMonoMass() nogil except +
        void setDiffAverageMass(double mass) nogil except +
        double getDiffAverageMass() nogil except +
        void setDiffMonoMass(double mass) nogil except +
        double getDiffMonoMass() nogil except +

        void setFormula(const String & composition) nogil except +
        String  getFormula() nogil except +
        void setDiffFormula(EmpiricalFormula & diff_formula) nogil except +
        EmpiricalFormula  getDiffFormula() nogil except +

        void setSynonyms(libcpp_set[ String ] & synonyms) nogil except +
        void addSynonym(const String & synonym) nogil except +
        libcpp_set[ String ] getSynonyms() nogil except +

        void setNeutralLossDiffFormulas(libcpp_vector[ EmpiricalFormula ] & diff_formulas) nogil except +
        libcpp_vector[ EmpiricalFormula ] getNeutralLossDiffFormulas() nogil except +

        void setNeutralLossMonoMasses(libcpp_vector[ double ] mono_masses) nogil except +
        libcpp_vector[ double ] getNeutralLossMonoMasses() nogil except +
        void setNeutralLossAverageMasses(libcpp_vector[ double ] average_masses) nogil except +
        libcpp_vector[ double ] getNeutralLossAverageMasses() nogil except +

        bool hasNeutralLoss() nogil except +
        bool isUserDefined() nogil except +

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
