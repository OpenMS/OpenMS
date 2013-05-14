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

        TargetedExperiment()                  nogil except +
        TargetedExperiment(TargetedExperiment)   nogil except +
        void clear(bool clear_meta_data)  nogil except +

        void sortTransitionsByProductMZ()  nogil except +

        TargetedExperiment operator+(TargetedExperiment)    nogil except +
        TargetedExperiment iadd(TargetedExperiment)   nogil except + # wrap-as:operator+=

        # cv list
        void setCVs(libcpp_vector[CV] cvs) nogil except +
        libcpp_vector[CV] getCVs() nogil except +
        void addCV(CV cv) nogil except +

        # contact list
        void setContacts(libcpp_vector[Contact] contacts) nogil except +
        libcpp_vector[Contact] getContacts() nogil except +
        void addContact(Contact contact) nogil except +

        # publication list
        void setPublications(libcpp_vector[Publication] publications) nogil except +
        libcpp_vector[Publication] getPublications() nogil except +
        void addPublication(Publication publication) nogil except +

        # target list
        void setTargetCVTerms(CVTermList cv_terms) nogil except +
        CVTermList getTargetCVTerms() nogil except +
        void addTargetCVTerm(CVTerm cv_term) nogil except +
        void setTargetMetaValue(String name, DataValue value) nogil except +

        # instrument list
        void setInstruments(libcpp_vector[TargetedExperiment_Instrument] instruments) nogil except +
        libcpp_vector[TargetedExperiment_Instrument] getInstruments() nogil except +
        void addInstrument(TargetedExperiment_Instrument instrument) nogil except +

        # software list
        void setSoftware(libcpp_vector[Software] software) nogil except +
        libcpp_vector[Software] getSoftware() nogil except +
        void addSoftware(Software software) nogil except +

        # protein list
        void setProteins(libcpp_vector[Protein] proteins) nogil except +
        libcpp_vector[Protein] getProteins() nogil except +
        Protein getProteinByRef(String ref) nogil except +
        void addProtein(Protein protein) nogil except +

        # compound list
        void setCompounds(libcpp_vector[Compound] rhs) nogil except +
        libcpp_vector[Compound] getCompounds() nogil except +
        void addCompound(Compound rhs) nogil except +

        void setPeptides(libcpp_vector[Peptide] rhs) nogil except +
        libcpp_vector[Peptide] getPeptides() nogil except +
        Peptide getPeptideByRef(String ref) nogil except +
        void addPeptide(Peptide rhs) nogil except +

        # set transition list
        void setTransitions(libcpp_vector[ReactionMonitoringTransition] transitions) nogil except +
        libcpp_vector[ReactionMonitoringTransition] getTransitions() nogil except +
        void addTransition(ReactionMonitoringTransition transition) nogil except +

        void setIncludeTargets(libcpp_vector[IncludeExcludeTarget] targets) nogil except +
        libcpp_vector[IncludeExcludeTarget] getIncludeTargets() nogil except +
        void addIncludeTarget(IncludeExcludeTarget target) nogil except +
        void setExcludeTargets(libcpp_vector[IncludeExcludeTarget] targets) nogil except +
        libcpp_vector[IncludeExcludeTarget] getExcludeTargets() nogil except +
        void addExcludeTarget(IncludeExcludeTarget target) nogil except +

        # sets the source files
        void setSourceFiles(libcpp_vector[SourceFile] source_files) nogil except +
        libcpp_vector[SourceFile] getSourceFiles() nogil except +
        void addSourceFile(SourceFile source_file) nogil except +

