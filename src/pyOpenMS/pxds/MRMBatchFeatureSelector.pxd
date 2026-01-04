from Types cimport *
from FeatureMap cimport *
from MRMFeatureSelector cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h>" namespace "OpenMS":

    cdef cppclass MRMBatchFeatureSelector:
        # wrap-doc:
        #   Batch processing wrapper for MRMFeatureSelector
        #
        #   This class enables iterative refinement of feature selection by applying
        #   MRMFeatureSelector multiple times with different parameter sets. Each iteration
        #   uses the output of the previous iteration as input, allowing for progressive
        #   filtering and optimization.
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
        #       # Configure selection parameters
        #       params1 = MRMFeatureSelector.SelectorParameters()
        #       params1.nn_threshold = 4
        #       params1.optimal_threshold = 0.5
        #       params1.score_weights = {"peak_apex_int": MRMFeatureSelector.LambdaScore.LINEAR}
        #
        #       # Second pass with stricter threshold
        #       params2 = MRMFeatureSelector.SelectorParameters()
        #       params2.optimal_threshold = 0.8
        #
        #       # Run batch selection
        #       selected = FeatureMap()
        #       MRMBatchFeatureSelector.batchMRMFeaturesScore(features, selected, [params1, params2])

        # Note: Constructor is deleted in C++, so we don't wrap it
        # All methods are static

# Static methods need to be declared separately
cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h>" namespace "OpenMS::MRMBatchFeatureSelector":

    void batchMRMFeatures(
        MRMFeatureSelector & feature_selector,
        FeatureMap & features,
        FeatureMap & selected_features,
        libcpp_vector[SelectorParameters] & parameters
    ) except + nogil # wrap-attach:MRMBatchFeatureSelector
        # wrap-doc:
        #   Apply feature selection iteratively with multiple parameter sets
        #
        #   Calls feature_selector.selectMRMFeature() repeatedly, using the result
        #   of each cycle as input for the next cycle.
        #
        #   :param feature_selector: MRMFeatureSelector instance to use
        #   :param features: Input FeatureMap
        #   :param selected_features: Output FeatureMap with selected features
        #   :param parameters: Vector of SelectorParameters for each iteration

    void batchMRMFeaturesScore(
        FeatureMap & features,
        FeatureMap & selected_features,
        libcpp_vector[SelectorParameters] & parameters
    ) except + nogil # wrap-attach:MRMBatchFeatureSelector
        # wrap-doc:
        #   Apply score-based feature selection iteratively
        #
        #   Convenience method that uses MRMFeatureSelectorScore internally.
        #
        #   :param features: Input FeatureMap
        #   :param selected_features: Output FeatureMap with selected features
        #   :param parameters: Vector of SelectorParameters for each iteration

    void batchMRMFeaturesQMIP(
        FeatureMap & features,
        FeatureMap & selected_features,
        libcpp_vector[SelectorParameters] & parameters
    ) except + nogil # wrap-attach:MRMBatchFeatureSelector
        # wrap-doc:
        #   Apply QMIP-based feature selection iteratively
        #
        #   Convenience method that uses MRMFeatureSelectorQMIP internally.
        #
        #   :param features: Input FeatureMap
        #   :param selected_features: Output FeatureMap with selected features
        #   :param parameters: Vector of SelectorParameters for each iteration
