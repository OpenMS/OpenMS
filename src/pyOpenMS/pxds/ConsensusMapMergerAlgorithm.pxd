from libcpp.map cimport map as libcpp_map
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from ConsensusMap cimport *
from ExperimentalDesign cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapMergerAlgorithm(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger
        ConsensusMapMergerAlgorithm() except + nogil
        ConsensusMapMergerAlgorithm(ConsensusMapMergerAlgorithm&) except + nogil
        void mergeProteinsAcrossFractionsAndReplicates(ConsensusMap& cmap,
                                                       const ExperimentalDesign& exp_design) except + nogil
        void mergeAllIDRuns(ConsensusMap& cmap) except + nogil
        void mergeProteinIDRuns(ConsensusMap& cmap,
                                const libcpp_map[unsigned int, unsigned int]& mapIdx_to_new_protIDRun) except + nogil
