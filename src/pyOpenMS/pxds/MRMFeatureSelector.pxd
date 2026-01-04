from Types cimport *
from String cimport *
from Feature cimport *
from FeatureMap cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>" namespace "OpenMS":

    # Note: VariableType and LambdaScore enums are not exposed due to
    # limitations with C++11 enum class in autowrap. Use integer values instead.

    cdef cppclass SelectorParameters "OpenMS::MRMFeatureSelector::SelectorParameters":
        # wrap-doc:
        #   Structure to configure MRMFeatureSelector parameters
        #
        #   Attributes:
        #   - nn_threshold: Number of nearest neighbors to include in optimization (default: 4)
        #   - locality_weight: Weight compounds by retention time proximity (default: False)
        #   - select_transition_group: Use component groups instead of components (default: True)
        #   - segment_window_length: Number of components in the sliding window (default: 8)
        #   - segment_step_length: Step size for sliding window (default: 4)
        #   - optimal_threshold: Cutoff for selection, 0 < x < 1 (default: 0.5)
        #
        #   Example:
        #
        #   .. code-block:: python
        #
        #       from pyopenms import *
        #
        #       params = SelectorParameters()
        #       params.nn_threshold = 6
        #       params.optimal_threshold = 0.7
        #       params.select_transition_group = False

        SelectorParameters()
        SelectorParameters(SelectorParameters &) # no-wrap

        Int nn_threshold
        bool locality_weight
        bool select_transition_group
        Int segment_window_length
        Int segment_step_length
        # Note: variable_type and score_weights not exposed due to enum class limitations
        double optimal_threshold


    # Note: MRMFeatureSelector base class is not wrapped directly because it is
    # an abstract class (has pure virtual optimize() method). Use the concrete
    # derived classes MRMFeatureSelectorQMIP or MRMFeatureSelectorScore instead.

    cdef cppclass MRMFeatureSelectorQMIP:
        # wrap-doc:
        #   Select MRM features using Quadratic Mixed Integer Programming (QMIP)
        #
        #   This class selects MRMFeatures based on relative retention time using a
        #   quadratic mixed integer programming formulation. The method optimizes
        #   feature selection while maintaining retention time consistency between
        #   neighboring transitions.
        #
        #   Example usage:
        #
        #   .. code-block:: python
        #
        #       from pyopenms import *
        #
        #       # Load features
        #       features = FeatureMap()
        #       FeatureXMLFile().load("features.featureXML", features)
        #
        #       # Configure parameters
        #       params = SelectorParameters()
        #       params.nn_threshold = 4
        #       params.optimal_threshold = 0.5
        #
        #       # Run QMIP selection
        #       selector = MRMFeatureSelectorQMIP()
        #       selected = FeatureMap()
        #       selector.selectMRMFeature(features, selected, params)

        MRMFeatureSelectorQMIP() except + nogil
        MRMFeatureSelectorQMIP(MRMFeatureSelectorQMIP &) except + nogil

        void selectMRMFeature(
            FeatureMap & features,
            FeatureMap & selected_filtered,
            SelectorParameters & parameters
        ) except + nogil
            # wrap-doc:
            #   Select MRM features using quadratic mixed integer programming
            #
            #   The features are sorted by retention time and split into segments with
            #   the given step and window length. The features are then selected based on
            #   relative retention time relationships using QMIP optimization.
            #
            #   :param features: Input FeatureMap
            #   :param selected_filtered: Output FeatureMap with selected features
            #   :param parameters: SelectorParameters configuration


    cdef cppclass MRMFeatureSelectorScore:
        # wrap-doc:
        #   Select MRM features using score-weighted linear programming
        #
        #   This class selects MRMFeatures based on a linear programming formulation
        #   where each possible transition is weighted by a user-defined score
        #   (typically retention time and peak intensity). Multiple scoring functions
        #   can be combined using the score_weights parameter.
        #
        #   Example usage:
        #
        #   .. code-block:: python
        #
        #       from pyopenms import *
        #
        #       # Load features
        #       features = FeatureMap()
        #       FeatureXMLFile().load("features.featureXML", features)
        #
        #       # Configure parameters
        #       params = SelectorParameters()
        #       params.nn_threshold = 4
        #       params.optimal_threshold = 0.5
        #
        #       # Run score-based selection
        #       selector = MRMFeatureSelectorScore()
        #       selected = FeatureMap()
        #       selector.selectMRMFeature(features, selected, params)

        MRMFeatureSelectorScore() except + nogil
        MRMFeatureSelectorScore(MRMFeatureSelectorScore &) except + nogil

        void selectMRMFeature(
            FeatureMap & features,
            FeatureMap & selected_filtered,
            SelectorParameters & parameters
        ) except + nogil
            # wrap-doc:
            #   Select MRM features using score-weighted linear programming
            #
            #   The features are sorted by retention time and split into segments with
            #   the given step and window length. The features are then selected based on
            #   score-weighted optimization.
            #
            #   :param features: Input FeatureMap
            #   :param selected_filtered: Output FeatureMap with selected features
            #   :param parameters: SelectorParameters configuration
