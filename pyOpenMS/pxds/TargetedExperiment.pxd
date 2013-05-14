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

