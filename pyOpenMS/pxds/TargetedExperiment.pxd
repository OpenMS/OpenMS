from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from ReactionMonitoringTransition cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>" namespace "OpenMS":

    cdef cppclass TargetedExperiment:

        TargetedExperiment()                  nogil except +
        TargetedExperiment(TargetedExperiment)   nogil except + #wrap-ignore
        void clear(bool clear_meta_data)

        libcpp_vector[ReactionMonitoringTransition] getTransitions()
        void setTransitions(libcpp_vector[ReactionMonitoringTransition] & transitions)
        void addTransition(ReactionMonitoringTransition & transition) 

        libcpp_vector[Peptide] getPeptides()
        Peptide getPeptideByRef(String & ref)
        void setPeptides(libcpp_vector[Peptide] & peptides)
        void addPeptide(Peptide & peptide)

        libcpp_vector[Protein] getProteins()
        Protein getProteinByRef(String & ref)
        void setProteins(libcpp_vector[Protein] & peptides)
        void addProtein(Protein & peptide)

        void sortTransitionsByProductMZ()


