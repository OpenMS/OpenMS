from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from CVTerm cimport *
from Residue cimport *
from Map cimport *
from DataValue cimport *
from CVTermList cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit":

    ctypedef enum RTUnit "OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit":
            # wrap-attach:
            #    RetentionTime
            SECOND,        # RT stored in seconds
            MINUTE,        # RT stored in minutes
            UNKNOWN,       # no stored annotation
            SIZE_OF_RTUNIT

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper::RetentionTime::RTType":

    ctypedef enum RTType "OpenMS::TargetedExperimentHelper::RetentionTime::RTType":
            # wrap-attach:
            #    RetentionTime
            LOCAL,            # undefined local chromatography
            NORMALIZED,       # standardized reference chromatography
            PREDICTED,        # predicted by referenced software
            HPINS,            # H-PINS "The de facto standard providing the retention times"
            IRT,              # iRT retention time standard
            UNKNOWN,          # no stored annotation
            SIZE_OF_RTTYPE

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper":

    cdef cppclass Configuration "OpenMS::TargetedExperimentHelper::Configuration":
        Configuration(Configuration) nogil except + #wrap-ignore
        String contact_ref
        String instrument_ref
        libcpp_vector[ CVTermList ] validations

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass CV:
        CV(CV) nogil except +
        CV(String new_id, String new_fullname, String new_version, String new_URI)  nogil except +

        String id
        String fullname
        String version
        String URI

    cdef cppclass Protein "OpenMS::TargetedExperimentHelper::Protein":
        Protein() nogil except +
        Protein(Protein) nogil except +

        String id
        String sequence

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass RetentionTime "OpenMS::TargetedExperimentHelper::RetentionTime":
        RetentionTime() nogil except +
        RetentionTime(RetentionTime) nogil except +

        String software_ref
        RTUnit retention_time_unit
        RTType retention_time_type

        bool isRTset() nogil except + 
        void setRT(double rt) nogil except +
        double getRT() nogil except +

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass Compound "OpenMS::TargetedExperimentHelper::Compound":
        Compound() nogil except +
        Compound(Compound) nogil except +
        bool operator==(Compound & rhs) nogil except +

        String id
        String molecular_formula
        String smiles_string
        double theoretical_mass
        libcpp_vector[ RetentionTime ] rts

        void setChargeState(int charge) nogil except +
        int getChargeState() nogil except +
        bool hasCharge() nogil except +
        double getRetentionTime() nogil except +
        bool hasRetentionTime() nogil except +
        RTType getRetentionTimeType() nogil except +
        RTUnit getRetentionTimeUnit() nogil except +

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass Peptide "OpenMS::TargetedExperimentHelper::Peptide":
        Peptide() nogil except +
        Peptide(Peptide) nogil except +

        # members
        libcpp_vector[RetentionTime] rts
        String id
        libcpp_vector[String] protein_refs
        CVTermList evidence
        String sequence
        libcpp_vector[TargetedExperiment_Modification] mods

        void setPeptideGroupLabel(String label) nogil except +
        String getPeptideGroupLabel() nogil except +

        void setChargeState(int charge) nogil except +
        int getChargeState() nogil except +
        bool hasCharge() nogil except +
        double getRetentionTime() nogil except +
        bool hasRetentionTime() nogil except +
        RTType getRetentionTimeType() nogil except +
        RTUnit getRetentionTimeUnit() nogil except +

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass Contact "OpenMS::TargetedExperimentHelper::Contact":
        Contact() nogil except +
        Contact(Contact) nogil except + #wrap-ignore
        String id
        bool operator==(Contact & rhs) nogil except +

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass Publication "OpenMS::TargetedExperimentHelper::Publication":
        Publication() nogil except +
        Publication(Publication) nogil except + #wrap-ignore
        String id
        bool operator==(Publication & rhs) nogil except +

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass TargetedExperiment_Instrument "OpenMS::TargetedExperimentHelper::Instrument":
        TargetedExperiment_Instrument() nogil except +
        TargetedExperiment_Instrument(TargetedExperiment_Instrument) nogil except + #wrap-ignore
        String id
        bool operator==(TargetedExperiment_Instrument & rhs) nogil except +

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass Prediction "OpenMS::TargetedExperimentHelper::Prediction":
        Prediction() nogil except +
        Prediction(Prediction) nogil except + #wrap-ignore
        bool operator==(Prediction & rhs) nogil except +

        String software_ref
        String contact_ref

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass TargetedExperiment_Interpretation "OpenMS::TargetedExperimentHelper::Interpretation":
        TargetedExperiment_Interpretation() nogil except +
        TargetedExperiment_Interpretation(TargetedExperiment_Interpretation) nogil except + #wrap-ignore

        unsigned char ordinal
        unsigned char rank
        ResidueType iontype

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass TraMLProduct "OpenMS::TargetedExperimentHelper::TraMLProduct":
        TraMLProduct() nogil except +
        TraMLProduct(TraMLProduct) nogil except + #wrap-ignore
        bool operator==(TraMLProduct & rhs) nogil except +

        void setMZ(double mz) nogil except +
        double getMZ() nogil except +

        void setChargeState(int charge) nogil except +
        int getChargeState() nogil except +
        bool hasCharge() nogil except +
        libcpp_vector[ Configuration ]  getConfigurationList() nogil except +
        void addConfiguration(Configuration configuration) nogil except +
        libcpp_vector[ TargetedExperiment_Interpretation ] getInterpretationList() nogil except +
        void addInterpretation(TargetedExperiment_Interpretation interpretation) nogil except +
        void resetInterpretations() nogil except +

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +
        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +
        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

# no support for nested classes yet in Cython
cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper::Peptide":

    cdef cppclass TargetedExperiment_Modification "OpenMS::TargetedExperimentHelper::Peptide::Modification":
        TargetedExperiment_Modification() nogil except +
        TargetedExperiment_Modification(TargetedExperiment_Modification) nogil except +

        # members
        double avg_mass_delta
        double mono_mass_delta
        int location
        int unimod_id

