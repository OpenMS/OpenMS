from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from ReactionMonitoringTransition cimport *
from TargetedExperimentHelper cimport *
from CVTermList cimport *
from IncludeExcludeTarget cimport *
from DataValue cimport *
from Software cimport *
from SourceFile cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>" namespace "OpenMS":

    cdef cppclass TargetedExperiment:

        TargetedExperiment() except + nogil  # TODO
        TargetedExperiment(TargetedExperiment &) except + nogil 

        bool operator==(TargetedExperiment) except + nogil 
        bool operator!=(TargetedExperiment) except + nogil 

        TargetedExperiment operator+(TargetedExperiment)    except + nogil 
        TargetedExperiment iadd(TargetedExperiment)   except + nogil  # wrap-as:operator+=

        void clear(bool clear_meta_data)  except + nogil 
        void sortTransitionsByProductMZ()  except + nogil 

        # cv list
        void setCVs(libcpp_vector[CV] cvs) except + nogil 
        libcpp_vector[CV] getCVs() except + nogil 
        void addCV(CV cv) except + nogil 

        # contact list
        void setContacts(libcpp_vector[Contact] contacts) except + nogil 
        libcpp_vector[Contact] getContacts() except + nogil 
        void addContact(Contact contact) except + nogil 

        # publication list
        void setPublications(libcpp_vector[Publication] publications) except + nogil 
        libcpp_vector[Publication] getPublications() except + nogil 
        void addPublication(Publication publication) except + nogil 

        # target list
        void setTargetCVTerms(CVTermList cv_terms) except + nogil 
        CVTermList getTargetCVTerms() except + nogil 
        void addTargetCVTerm(CVTerm cv_term) except + nogil 
        void setTargetMetaValue(String name, DataValue value) except + nogil 

        # instrument list
        void setInstruments(libcpp_vector[TargetedExperiment_Instrument] instruments) except + nogil 
        libcpp_vector[TargetedExperiment_Instrument] getInstruments() except + nogil 
        void addInstrument(TargetedExperiment_Instrument instrument) except + nogil 

        # software list
        void setSoftware(libcpp_vector[Software] software) except + nogil 
        libcpp_vector[Software] getSoftware() except + nogil 
        void addSoftware(Software software) except + nogil 

        # protein list
        void setProteins(libcpp_vector[Protein] proteins) except + nogil 
        libcpp_vector[Protein] getProteins() except + nogil 
        Protein getProteinByRef(String ref) except + nogil 
        bool hasProtein(String ref) except + nogil 
        void addProtein(Protein protein) except + nogil 

        # compound list
        void setCompounds(libcpp_vector[Compound] rhs) except + nogil 
        libcpp_vector[Compound] getCompounds() except + nogil 
        void addCompound(Compound rhs) except + nogil 
        bool hasCompound(String ref) except + nogil 
        Compound getCompoundByRef(String ref) except + nogil 

        # peptide list
        void setPeptides(libcpp_vector[Peptide] rhs) except + nogil 
        libcpp_vector[Peptide] getPeptides() except + nogil 
        bool hasPeptide(String ref) except + nogil 
        Peptide getPeptideByRef(String ref) except + nogil 
        void addPeptide(Peptide rhs) except + nogil 

        # set transition list
        void setTransitions(libcpp_vector[ReactionMonitoringTransition] transitions) except + nogil 
        libcpp_vector[ReactionMonitoringTransition] getTransitions() except + nogil 
        void addTransition(ReactionMonitoringTransition transition) except + nogil 

        void setIncludeTargets(libcpp_vector[IncludeExcludeTarget] targets) except + nogil 
        libcpp_vector[IncludeExcludeTarget] getIncludeTargets() except + nogil 
        void addIncludeTarget(IncludeExcludeTarget target) except + nogil 
        void setExcludeTargets(libcpp_vector[IncludeExcludeTarget] targets) except + nogil 
        libcpp_vector[IncludeExcludeTarget] getExcludeTargets() except + nogil 
        void addExcludeTarget(IncludeExcludeTarget target) except + nogil 

        # sets the source files
        void setSourceFiles(libcpp_vector[SourceFile] source_files) except + nogil 
        libcpp_vector[SourceFile] getSourceFiles() except + nogil 
        void addSourceFile(SourceFile source_file) except + nogil 

        bool containsInvalidReferences() except + nogil 
