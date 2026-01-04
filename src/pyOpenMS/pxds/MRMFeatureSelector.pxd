from Types cimport *
from String cimport *
from Feature cimport *
from FeatureMap cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>" namespace "OpenMS::MRMFeatureSelector":

    cdef enum VariableType "OpenMS::MRMFeatureSelector::VariableType":
        # wrap-attach:
        #   MRMFeatureSelector
        INTEGER
        CONTINUOUS

    cdef enum LambdaScore "OpenMS::MRMFeatureSelector::LambdaScore":
        # wrap-attach:
        #   MRMFeatureSelector
        LINEAR
        INVERSE
        LOG
        INVERSE_LOG
        INVERSE_LOG10

    cdef cppclass SelectorParameters "OpenMS::MRMFeatureSelector::SelectorParameters":
        # wrap-attach:
        #   MRMFeatureSelector
        # wrap-doc:
        #   Structure to configure MRMFeatureSelector parameters
        #
        #   Attributes:
        #   - nn_threshold: Number of nearest neighbors to include in optimization
        #   - locality_weight: Weight compounds by retention time proximity
        #   - select_transition_group: Use component groups instead of components
        #   - segment_window_length: Number of components in the sliding window
        #   - segment_step_length: Step size for sliding window
        #   - variable_type: INTEGER for exact solution, CONTINUOUS for relaxed
        #   - optimal_threshold: Cutoff for selection (0 < x < 1)
        #   - score_weights: Map of score names to lambda transformations

        SelectorParameters() except + nogil
        SelectorParameters(SelectorParameters &) except + nogil

        Int nn_threshold
        bool locality_weight
        bool select_transition_group
        Int segment_window_length
        Int segment_step_length
        VariableType variable_type
        double optimal_threshold
        libcpp_map[String, LambdaScore] score_weights


cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureSelector:
        # wrap-doc:
        #   Base class for selection of MRM Features through Linear Programming optimization
        #
        #   This class provides the framework for optimal feature selection in MRM/SRM experiments
        #   using Linear Programming (LP). The key idea is to select the best peak for each transition
        #   while maintaining consistency with neighboring transitions based on retention time relationships.
        #
        #   Two derived implementations are provided:
        #   - MRMFeatureSelectorQMIP: Uses Quadratic Mixed Integer Programming based on relative RT
        #   - MRMFeatureSelectorScore: Uses score-weighted linear programming
        #
        #   Algorithm:
        #   1. Features are sorted by retention time
        #   2. The RT range is divided into overlapping segments (sliding window)
        #   3. For each segment, an LP problem is formulated and solved
        #   4. Solutions from segments are merged to produce final selection

        MRMFeatureSelector() except + nogil
        MRMFeatureSelector(MRMFeatureSelector &) except + nogil

        void selectMRMFeature(
            FeatureMap & features,
            FeatureMap & selected_filtered,
            SelectorParameters & parameters
        ) except + nogil
            # wrap-doc:
            #   Select MRM features using linear programming optimization
            #
            #   The features are sorted by retention time and split into segments with
            #   the given step and window length. The features are then selected based on
            #   the results of the optimize() method applied to each segment.
            #
            #   :param features: Input FeatureMap
            #   :param selected_filtered: Output FeatureMap with selected features
            #   :param parameters: SelectorParameters configuration


    cdef cppclass MRMFeatureSelectorQMIP(MRMFeatureSelector):
        # wrap-inherits:
        #   MRMFeatureSelector
        # wrap-doc:
        #   Select MRM features using Quadratic Mixed Integer Programming (QMIP)
        #
        #   This class selects MRMFeatures based on relative retention time using a
        #   quadratic mixed integer programming formulation. The method optimizes
        #   feature selection while maintaining retention time consistency between
        #   neighboring transitions.

        MRMFeatureSelectorQMIP() except + nogil
        MRMFeatureSelectorQMIP(MRMFeatureSelectorQMIP &) except + nogil


    cdef cppclass MRMFeatureSelectorScore(MRMFeatureSelector):
        # wrap-inherits:
        #   MRMFeatureSelector
        # wrap-doc:
        #   Select MRM features using score-weighted linear programming
        #
        #   This class selects MRMFeatures based on a linear programming formulation
        #   where each possible transition is weighted by a user-defined score
        #   (typically retention time and peak intensity). Multiple scoring functions
        #   can be combined using the score_weights parameter.
        #
        #   Supported lambda transformations for scores:
        #   - LINEAR: Use score as-is
        #   - INVERSE: 1/score
        #   - LOG: log(score)
        #   - INVERSE_LOG: 1/log(score)
        #   - INVERSE_LOG10: 1/log10(score)

        MRMFeatureSelectorScore() except + nogil
        MRMFeatureSelectorScore(MRMFeatureSelectorScore &) except + nogil
