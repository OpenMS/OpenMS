from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS":
    
    cdef cppclass ResidueModification "OpenMS::ResidueModification":
        ResidueModification() nogil except +
        ResidueModification(ResidueModification) nogil except +
        void setId(String & id_) nogil except +
        String  getId() nogil except +
        void setFullId(String & full_id) nogil except +
        String  getFullId() nogil except +
        void setUniModAccession(String & id_) nogil except +
        String  getUniModAccession() nogil except +
        void setPSIMODAccession(String & id_) nogil except +
        String  getPSIMODAccession() nogil except +
        void setFullName(String & full_name) nogil except +
        String  getFullName() nogil except +
        void setName(String & name) nogil except +
        String  getName() nogil except +
        void setTermSpecificity(Term_Specificity term_spec) nogil except +
        void setTermSpecificity(String & name) nogil except +
        Term_Specificity getTermSpecificity() nogil except +
        String getTermSpecificityName(Term_Specificity ) nogil except +
        void setOrigin(String & origin) nogil except +
        String  getOrigin() nogil except +
        void setSourceClassification(String & classification) nogil except +
        void setSourceClassification(Source_Classification classification) nogil except +
        Source_Classification getSourceClassification() nogil except +
        String getSourceClassificationName(Source_Classification classification) nogil except +
        void setAverageMass(DoubleReal mass) nogil except +
        DoubleReal getAverageMass() nogil except +
        void setMonoMass(DoubleReal mass) nogil except +
        DoubleReal getMonoMass() nogil except +
        void setDiffAverageMass(DoubleReal mass) nogil except +
        DoubleReal getDiffAverageMass() nogil except +
        void setDiffMonoMass(DoubleReal mass) nogil except +
        DoubleReal getDiffMonoMass() nogil except +
        void setFormula(String & composition) nogil except +
        String  getFormula() nogil except +
        void setDiffFormula(EmpiricalFormula & diff_formula) nogil except +
        EmpiricalFormula  getDiffFormula() nogil except +
        void setSynonyms(libcpp_set[ String ] & synonyms) nogil except +
        void addSynonym(String & synonym) nogil except +
        libcpp_set[ String ]  getSynonyms() nogil except +
        void setNeutralLossDiffFormula(EmpiricalFormula & loss) nogil except +
        EmpiricalFormula  getNeutralLossDiffFormula() nogil except +
        void setNeutralLossMonoMass(DoubleReal mono_mass) nogil except +
        DoubleReal getNeutralLossMonoMass() nogil except +
        void setNeutralLossAverageMass(DoubleReal average_mass) nogil except +
        DoubleReal getNeutralLossAverageMass() nogil except +
        bool hasNeutralLoss() nogil except +
        bool operator==(ResidueModification & modification) nogil except +
        bool operator!=(ResidueModification & modification) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS::ResidueModification":
    cdef enum Term_Specificity "OpenMS::ResidueModification::Term_Specificity":
        #wrap-attach:
        #    ResidueModification
        ANYWHERE
        C_TERM
        N_TERM
        NUMBER_OF_TERM_SPECIFICITY

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS::ResidueModification":
    cdef enum Source_Classification "OpenMS::ResidueModification::Source_Classification":
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

