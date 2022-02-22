from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from MSExperiment cimport *
from LPWrapper cimport *
from String cimport *
from PrecursorIonSelectionPreprocessing cimport *
from PSProteinInference cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/PSLPFormulation.h>" namespace "OpenMS":
    
    cdef cppclass PSLPFormulation(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        PSLPFormulation() nogil except +
        PSLPFormulation(PSLPFormulation &) nogil except + # compiler

        void createAndSolveILPForKnownLCMSMapFeatureBased(
            FeatureMap & features,
            MSExperiment & experiment,
            libcpp_vector[ IndexTriple ] & variable_indices,
            libcpp_vector[ libcpp_vector[ libcpp_pair[ size_t, size_t ] ] ] & mass_ranges, 
            libcpp_set[ int ] & charges_set,
            UInt ms2_spectra_per_rt_bin,
            libcpp_vector[ int ] & solution_indices) nogil except + 
                # wrap-doc:
                #   Encode ILP formulation for a given LC-MS map, but unknown protein sample
                #   -----
                #   :param features: FeatureMap with all possible precursors
                #   :param experiment: Input raw data
                #   :param variable_indices: Assignment of feature indices and ILP variables
                #   :param mass_ranges: Feature borders as indices in the raw data
                #   :param charges_set: Allowed charge states
                #   :param ms2_spectra_per_rt_bin: Allowed number of precursors per rt bin
                #   :param solution_indices: Indices of ILP variables that are in the optimal solution

        void createAndSolveILPForInclusionListCreation(
            PrecursorIonSelectionPreprocessing & preprocessing,
            UInt ms2_spectra_per_rt_bin, UInt max_list_size, FeatureMap & precursors, bool solve_ILP) nogil except + # wrap-doc:Find a set of precursors, so that the protein coverage is maximal and that the number of precursors per bin is not exceeded

        void createAndSolveCombinedLPForKnownLCMSMapFeatureBased(FeatureMap & features, MSExperiment & experiment, 
            libcpp_vector[ IndexTriple ] & variable_indices, libcpp_vector[ int ] & solution_indices, 
            libcpp_vector[ libcpp_vector[ libcpp_pair[ size_t, size_t ] ] ] & mass_ranges, 
            libcpp_set[ Int ] & charges_set, UInt ms2_spectra_per_rt_bin, Size step_size, bool sequential_order) nogil except + # wrap-ignore

        void updateStepSizeConstraint(Size iteration, UInt step_size) nogil except +

        void updateFeatureILPVariables(
          FeatureMap & new_features,
          libcpp_vector[ IndexTriple ] & variable_indices,
          libcpp_map[ Size, libcpp_vector[ String ] ] & feature_constraints_map) nogil except + # wrap-ignore

        void updateRTConstraintsForSequentialILP(Size & rt_index, UInt ms2_spectra_per_rt_bin, Size max_rt_index) nogil except +

        void updateCombinedILP(FeatureMap & features,
                               PrecursorIonSelectionPreprocessing & preprocessed_db,
                               libcpp_vector[ IndexTriple ] & variable_indices,
                               libcpp_vector[ String ] & new_protein_accs,
                               libcpp_vector[ String ] & protein_accs,
                               PSProteinInference & prot_inference,
                               Size & variable_counter,
                               libcpp_map[ String, libcpp_vector[ size_t ] ] & protein_feature_map,
                               Feature & new_feature, libcpp_map[ String, Size ] & protein_variable_index_map,
                               libcpp_map[ String, libcpp_set[ String ] ] & prot_id_counter) nogil except + # wrap-ignore

        void solveILP(libcpp_vector[ int ] & solution_indices) nogil except + # wrap-doc:Solve the ILP

        void setLPSolver(SOLVER solver) nogil except +

        SOLVER getLPSolver() nogil except +

cdef extern from "<OpenMS/ANALYSIS/TARGETED/PSLPFormulation.h>" namespace "OpenMS::PSLPFormulation":
    
    cdef cppclass IndexTriple "OpenMS::PSLPFormulation::IndexTriple":
        IndexTriple() nogil except +
        IndexTriple(IndexTriple) nogil except + #wrap-ignore
        Size feature
        Int scan
        Size variable
        double rt_probability
        double signal_weight
        String prot_acc

