from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from ReactionMonitoringTransition cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>" namespace "OpenMS":

    cdef cppclass TargetedExperiment:

        TargetedExperiment()                  nogil except +
        TargetedExperiment(TargetedExperiment)   nogil except +
        void clear(bool clear_meta_data)  nogil except +

        libcpp_vector[ReactionMonitoringTransition] getTransitions()  nogil except +
        void setTransitions(libcpp_vector[ReactionMonitoringTransition] transitions)  nogil except +
        void addTransition(ReactionMonitoringTransition transition)   nogil except +

        libcpp_vector[Peptide] getPeptides()  nogil except +
        Peptide getPeptideByRef(String ref)  nogil except +
        void setPeptides(libcpp_vector[Peptide] peptides)  nogil except +
        void addPeptide(Peptide peptide)  nogil except +

        libcpp_vector[Protein] getProteins()  nogil except +
        Protein getProteinByRef(String ref)  nogil except +
        void setProteins(libcpp_vector[Protein] peptides)  nogil except +
        void addProtein(Protein peptide)  nogil except +

        void sortTransitionsByProductMZ()  nogil except +

        TargetedExperiment operator+(TargetedExperiment)    nogil except +
        TargetedExperiment iadd(TargetedExperiment)   nogil except + # wrap-as:operator+=

"""
        # cv list
        void setCVs(libcpp_vector[CV] & cvs) nogil except +

        libcpp_vector[CV] & getCVs() nogil except +

        void addCV(CV & cv) nogil except +
        # contact list
        void setContacts(libcpp_vector[Contact] & contacts) nogil except +

        libcpp_vector[Contact] & getContacts() nogil except +

        void addContact(Contact & contact) nogil except +

        # publication list
        void setPublications(libcpp_vector[Publication] & publications) nogil except +

        libcpp_vector[Publication] & getPublications() nogil except +

        void addPublication(Publication & publication) nogil except +

        # target list
        void setTargetCVTerms(CVTermList & cv_terms) nogil except +

        CVTermList & getTargetCVTerms() nogil except +

        void addTargetCVTerm(CVTerm & cv_term) nogil except +

        void setTargetMetaValue(String & name, DataValue & value) nogil except +

        # instrument list
        void setInstruments(libcpp_vector[Instrument] & instruments) nogil except +

        libcpp_vector[Instrument] & getInstruments() nogil except +

        void addInstrument(Instrument & instrument) nogil except +

        # software list
        void setSoftware(libcpp_vector[Software] & software) nogil except +

        libcpp_vector[Software] & getSoftware() nogil except +

        void addSoftware(Software & software) nogil except +

        # protein list
        void setProteins(libcpp_vector[Protein] & proteins) nogil except +

        libcpp_vector[Protein] & getProteins() nogil except +

        Protein & getProteinByRef(String & ref) nogil except +

        void addProtein(Protein & protein) nogil except +

        # compound list
        void setCompounds(libcpp_vector[Compound] & rhs) nogil except +

        libcpp_vector[Compound] & getCompounds() nogil except +

        void addCompound(Compound & rhs) nogil except +

        void setPeptides(libcpp_vector[Peptide] & rhs) nogil except +

        libcpp_vector[Peptide] & getPeptides() nogil except +

        Peptide & getPeptideByRef(String & ref) nogil except +

        void addPeptide(Peptide & rhs) nogil except +

        # set transition list
        void setTransitions(libcpp_vector[ReactionMonitoringTransition] & transitions) nogil except +

        # returns the transition list
        libcpp_vector[ReactionMonitoringTransition] & getTransitions() nogil except +

        # adds a transition to the list
        void addTransition(ReactionMonitoringTransition & transition) nogil except +

        void setIncludeTargets(libcpp_vector[IncludeExcludeTarget] & targets) nogil except +

        libcpp_vector[IncludeExcludeTarget] & getIncludeTargets() nogil except +

        void addIncludeTarget(IncludeExcludeTarget & target) nogil except +

        void setExcludeTargets(libcpp_vector[IncludeExcludeTarget] & targets) nogil except +

        libcpp_vector[IncludeExcludeTarget] & getExcludeTargets() nogil except +

        void addExcludeTarget(IncludeExcludeTarget & target) nogil except +

        # sets the source files
        void setSourceFiles(libcpp_vector[SourceFile] & source_files) nogil except +

        # returns the source file list
        libcpp_vector[SourceFile] & getSourceFiles() nogil except +

        # adds a source file to the list
        void addSourceFile(SourceFile & source_file) nogil except +
        #@}

        #@name Sorting peaks
        #@{
        /**
          @brief Lexicographically sorts the transitions by their product m/z.
        */
        void sortTransitionsByProductMZ() nogil except +
        #@}

    protected:

        void createProteinReferenceMap_() nogil except +

        void createPeptideReferenceMap_() nogil except +

        libcpp_vector[CV] cvs_ nogil except +

        libcpp_vector[Contact] contacts_ nogil except +

        libcpp_vector[Publication] publications_ nogil except +

        libcpp_vector[Instrument] instruments_ nogil except +

        CVTermList targets_ nogil except +

        libcpp_vector[Software] software_ nogil except +

        libcpp_vector[Protein] proteins_ nogil except +

        libcpp_vector[Compound] compounds_ nogil except +

        libcpp_vector[Peptide] peptides_ nogil except +

        libcpp_vector[ReactionMonitoringTransition] transitions_ nogil except +

        libcpp_vector[IncludeExcludeTarget] include_targets_ nogil except +

        libcpp_vector[IncludeExcludeTarget] exclude_targets_ nogil except +

        libcpp_vector[SourceFile] source_files_ nogil except +

        ProteinReferenceMapType protein_reference_map_ nogil except +

        bool protein_reference_map_dirty_ nogil except +

        PeptideReferenceMapType peptide_reference_map_ nogil except +

        bool peptide_reference_map_dirty_ nogil except +

      } nogil except +


      namespace TargetedExperimentHelper
      {
      } # namespace TargetedExperimentHelper


    } # namespace OpenMS






    TODO: Found function in cpp but not in pxd: function public setCVs
    TODO: Found function in cpp but not in pxd: function public getCVs
    TODO: Found function in cpp but not in pxd: function public addCV
    TODO: Found function in cpp but not in pxd: function public setContacts
    TODO: Found function in cpp but not in pxd: function public getContacts
    TODO: Found function in cpp but not in pxd: function public addContact
    TODO: Found function in cpp but not in pxd: function public setPublications
    TODO: Found function in cpp but not in pxd: function public getPublications
    TODO: Found function in cpp but not in pxd: function public addPublication
    TODO: Found function in cpp but not in pxd: function public setTargetCVTerms
    TODO: Found function in cpp but not in pxd: function public getTargetCVTerms
    TODO: Found function in cpp but not in pxd: function public addTargetCVTerm
    TODO: Found function in cpp but not in pxd: function public setTargetMetaValue
    TODO: Found function in cpp but not in pxd: function public setInstruments
    TODO: Found function in cpp but not in pxd: function public getInstruments
    TODO: Found function in cpp but not in pxd: function public addInstrument
    TODO: Found function in cpp but not in pxd: function public setSoftware
    TODO: Found function in cpp but not in pxd: function public getSoftware
    TODO: Found function in cpp but not in pxd: function public addSoftware
    TODO: Found function in cpp but not in pxd: function public setCompounds
    TODO: Found function in cpp but not in pxd: function public getCompounds
    TODO: Found function in cpp but not in pxd: function public addCompound
    TODO: Found function in cpp but not in pxd: function public setIncludeTargets
    TODO: Found function in cpp but not in pxd: function public getIncludeTargets
    TODO: Found function in cpp but not in pxd: function public addIncludeTarget
    TODO: Found function in cpp but not in pxd: function public setExcludeTargets
    TODO: Found function in cpp but not in pxd: function public getExcludeTargets
    TODO: Found function in cpp but not in pxd: function public addExcludeTarget
    TODO: Found function in cpp but not in pxd: function public setSourceFiles
    TODO: Found function in cpp but not in pxd: function public getSourceFiles
    TODO: Found function in cpp but not in pxd: function public addSourceFile
    """
