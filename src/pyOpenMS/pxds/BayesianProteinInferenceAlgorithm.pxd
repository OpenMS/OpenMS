from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ExperimentalDesign cimport *
from PeptideIdentification cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>" namespace "OpenMS":

    cdef cppclass BayesianProteinInferenceAlgorithm(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        #
        # wrap-doc:
        #     Performs a Bayesian protein inference on Protein/Peptide identifications or ConsensusMap.
        #     - Filters for best n PSMs per spectrum.
        #     - Calculates and filters for best peptide per spectrum.
        #     - Builds a k-partite graph from the structures.
        #     - Finds and splits into connected components by DFS
        #     - Extends the graph by adding layers from indist. protein groups, peptides with the same parents and optionally
        #       some additional layers (peptide sequence, charge, replicate -> extended model = experimental)
        #     - Builds a factor graph representation of a Bayesian network using the Evergreen library
        #       See model param section. It is based on the Fido noisy-OR model with an option for
        #       regularizing the number of proteins per peptide.
        #     - Performs loopy belief propagation on the graph and queries protein, protein group and/or peptide posteriors
        #       See loopy_belief_propagation param section.
        #     - Learns best parameters via grid search if the parameters were not given in the param section.
        #     - Writes posteriors to peptides and/or proteins and adds indistinguishable protein groups to the underlying
        #       data structures.
        #     - Can make use of OpenMP to parallelize over connected components.

        BayesianProteinInferenceAlgorithm() nogil except +

        BayesianProteinInferenceAlgorithm(BayesianProteinInferenceAlgorithm) nogil except + #wrap-ignore

        BayesianProteinInferenceAlgorithm(unsigned int debug_lvl) nogil except +


        void inferPosteriorProbabilities(libcpp_vector[ ProteinIdentification ] & proteinIDs, 
                                         libcpp_vector[ PeptideIdentification ] & peptideIDs) nogil except +

        void inferPosteriorProbabilities(libcpp_vector[ ProteinIdentification ] & proteinIDs, 
                                         libcpp_vector[ PeptideIdentification ] & peptideIDs, 
                                         ExperimentalDesign exp_des) nogil except +

        # not sure if that can be wrapped at all? 
        # void inferPosteriorProbabilities(ConsensusMap & cmap, bool greedy_group_resolution, boost::optional[ ExperimentalDesign ] exp_des) nogil except +
