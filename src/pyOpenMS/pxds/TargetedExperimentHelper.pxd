from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from CVTerm cimport *
from Residue cimport *
from DataValue cimport *
from CVTermList cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit":

    ctypedef enum RTUnit "OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit":
            # wrap-attach:
            #   RetentionTime
            SECOND,        # RT stored in seconds
            MINUTE,        # RT stored in minutes
            UNKNOWN,       # no stored annotation
            SIZE_OF_RTUNIT

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper::RetentionTime::RTType":

    ctypedef enum RTType "OpenMS::TargetedExperimentHelper::RetentionTime::RTType":
            # wrap-attach:
            #   RetentionTime
            LOCAL,            # undefined local chromatography
            NORMALIZED,       # standardized reference chromatography
            PREDICTED,        # predicted by referenced software
            HPINS,            # H-PINS "The de facto standard providing the retention times"
            IRT,              # iRT retention time standard
            UNKNOWN,          # no stored annotation
            SIZE_OF_RTTYPE

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper":

    cdef cppclass Configuration(CVTermList):
        # wrap-inherits:
        #   CVTermList

        Configuration(Configuration) except + nogil  #wrap-ignore
        String contact_ref
        String instrument_ref
        libcpp_vector[ CVTermList ] validations


    cdef cppclass CV:
        CV(CV &) except + nogil 
        CV(String new_id, String new_fullname, String new_version, String new_URI)  except + nogil 

        String id
        String fullname
        String version
        String URI

    cdef cppclass Protein(CVTermList):
        # wrap-inherits:
        #   CVTermList

        Protein() except + nogil 
        Protein(Protein &) except + nogil 

        String id
        String sequence


    cdef cppclass RetentionTime(CVTermList):
        # wrap-inherits:
        #   CVTermList

        RetentionTime() except + nogil 
        RetentionTime(RetentionTime &) except + nogil 

        String software_ref
        RTUnit retention_time_unit
        RTType retention_time_type

        bool isRTset() except + nogil  
        void setRT(double rt) except + nogil 
        double getRT() except + nogil 

    cdef cppclass Compound(CVTermList):
        # wrap-inherits:
        #   CVTermList

        Compound() except + nogil 
        Compound(Compound &) except + nogil 
        bool operator==(Compound & rhs) except + nogil 

        String id
        String molecular_formula
        String smiles_string
        double theoretical_mass
        libcpp_vector[ RetentionTime ] rts

        void setChargeState(int charge) except + nogil  # wrap-doc:Sets the peptide or compound charge state
        int getChargeState() except + nogil  # wrap-doc:Returns the peptide or compound charge state
        bool hasCharge() except + nogil  # wrap-doc:Whether peptide or compound has set charge state
        double getRetentionTime() except + nogil  # wrap-doc:Gets compound or peptide retention time
        bool hasRetentionTime() except + nogil  # wrap-doc:Check whether compound or peptide has an annotated retention time
        RTType getRetentionTimeType() except + nogil  # wrap-doc:Get compound or peptide retentiontime type
        RTUnit getRetentionTimeUnit() except + nogil  # wrap-doc:Get compound or peptide retentiontime type

    cdef cppclass Peptide(CVTermList):
        # wrap-inherits:
        #   CVTermList

        Peptide() except + nogil 
        Peptide(Peptide &) except + nogil 

        # members
        libcpp_vector[RetentionTime] rts
        String id
        libcpp_vector[String] protein_refs
        CVTermList evidence
        String sequence
        libcpp_vector[TargetedExperiment_Modification] mods

        void setPeptideGroupLabel(String label) except + nogil  # wrap-doc:Sets the peptide group label
        String getPeptideGroupLabel() except + nogil  # wrap-doc:Get the peptide group label

        void setChargeState(int charge) except + nogil  # wrap-doc:Sets the peptide or compound charge states
        int getChargeState() except + nogil  # wrap-doc:Returns the peptide or compound charge state
        bool hasCharge() except + nogil  # wrap-doc:Whether product has set charge state
        double getRetentionTime() except + nogil  # wrap-doc:Gets compound or peptide retention time
        bool hasRetentionTime() except + nogil  # wrap-doc:Gets compound or peptide retention time
        RTType getRetentionTimeType() except + nogil  # wrap-doc:Get compound or peptide retentiontime type
        RTUnit getRetentionTimeUnit() except + nogil  # wrap-doc:Get compound or peptide retentiontime unit (minute/seconds)

    cdef cppclass Contact(CVTermList):
        # wrap-inherits:
        #   CVTermList

        Contact() except + nogil 
        Contact(Contact &) except + nogil  # compiler
        String id
        bool operator==(Contact & rhs) except + nogil 

    cdef cppclass Publication(CVTermList):
        # wrap-inherits:
        #   CVTermList

        Publication() except + nogil 
        Publication(Publication &) except + nogil  # compiler
        String id
        bool operator==(Publication & rhs) except + nogil 


    cdef cppclass TargetedExperiment_Instrument "OpenMS::TargetedExperimentHelper::Instrument":
        TargetedExperiment_Instrument() except + nogil 
        TargetedExperiment_Instrument(TargetedExperiment_Instrument &) except + nogil  # compiler
        String id
        bool operator==(TargetedExperiment_Instrument & rhs) except + nogil 

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  except + nogil 
        void replaceCVTerm(CVTerm & term)               except + nogil 
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) except + nogil 
        void replaceCVTerms(libcpp_map[String, libcpp_vector[CVTerm] ] cv_term_map) except + nogil 
        void consumeCVTerms(libcpp_map[String, libcpp_vector[CVTerm] ] cv_term_map) except + nogil 
        libcpp_map[String, libcpp_vector[CVTerm] ] getCVTerms() except + nogil 
        void addCVTerm(CVTerm & term)                   except + nogil 
        bool hasCVTerm(String accession)  except + nogil 
        bool empty()                      except + nogil 

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) except + nogil 
        void getKeys(libcpp_vector[unsigned int] & keys) except + nogil  # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) except + nogil 
        DataValue getMetaValue(String) except + nogil 
        void setMetaValue(unsigned int, DataValue) except + nogil 
        void setMetaValue(String, DataValue) except + nogil 
        bool metaValueExists(String) except + nogil 
        bool metaValueExists(unsigned int) except + nogil 
        void removeMetaValue(String) except + nogil 
        void removeMetaValue(unsigned int) except + nogil 

    cdef cppclass Prediction(CVTermList):
        # wrap-inherits:
        #   CVTermList

        Prediction() except + nogil 
        Prediction(Prediction &) except + nogil  # compiler
        bool operator==(Prediction & rhs) except + nogil 

        String software_ref
        String contact_ref


    cdef cppclass TargetedExperiment_Interpretation "OpenMS::TargetedExperimentHelper::Interpretation":
        TargetedExperiment_Interpretation() except + nogil 
        TargetedExperiment_Interpretation(TargetedExperiment_Interpretation &) except + nogil  # compiler

        unsigned char ordinal
        unsigned char rank
        ResidueType iontype

        # CVTermList:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  except + nogil 
        void replaceCVTerm(CVTerm & term)               except + nogil 
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) except + nogil 
        void replaceCVTerms(libcpp_map[String, libcpp_vector[CVTerm] ] cv_term_map) except + nogil 
        void consumeCVTerms(libcpp_map[String, libcpp_vector[CVTerm] ] cv_term_map) except + nogil 
        libcpp_map[String, libcpp_vector[CVTerm] ] getCVTerms() except + nogil 
        void addCVTerm(CVTerm & term)                   except + nogil 
        bool hasCVTerm(String accession)  except + nogil 
        bool empty()                      except + nogil 

        # MetaInfoInterface:
        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) except + nogil 
        void getKeys(libcpp_vector[unsigned int] & keys) except + nogil  # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) except + nogil 
        DataValue getMetaValue(String) except + nogil 
        void setMetaValue(unsigned int, DataValue) except + nogil 
        void setMetaValue(String, DataValue) except + nogil 
        bool metaValueExists(String) except + nogil 
        bool metaValueExists(unsigned int) except + nogil 
        void removeMetaValue(String) except + nogil 
        void removeMetaValue(unsigned int) except + nogil 

    cdef cppclass TraMLProduct(CVTermList):
        # wrap-inherits:
        #   CVTermList

        TraMLProduct() except + nogil 
        TraMLProduct(TraMLProduct &) except + nogil  # compiler
        bool operator==(TraMLProduct & rhs) except + nogil 

        void setMZ(double mz) except + nogil 
        double getMZ() except + nogil 

        void setChargeState(int charge) except + nogil 
        int getChargeState() except + nogil 
        bool hasCharge() except + nogil 
        libcpp_vector[ Configuration ]  getConfigurationList() except + nogil 
        void addConfiguration(Configuration configuration) except + nogil 
        libcpp_vector[ TargetedExperiment_Interpretation ] getInterpretationList() except + nogil 
        void addInterpretation(TargetedExperiment_Interpretation interpretation) except + nogil 
        void resetInterpretations() except + nogil 


# no support for nested classes yet in Cython
cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper::Peptide":

    cdef cppclass TargetedExperiment_Modification "OpenMS::TargetedExperimentHelper::Peptide::Modification":
        TargetedExperiment_Modification() except + nogil 
        TargetedExperiment_Modification(TargetedExperiment_Modification &) except + nogil 

        # members
        double avg_mass_delta
        double mono_mass_delta
        int location
        int unimod_id

