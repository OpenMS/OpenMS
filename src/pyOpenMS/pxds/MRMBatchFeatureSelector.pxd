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
        #   Note: The static batch methods are not directly wrapped due to limitations
        #   with nested type handling. For iterative selection, use MRMFeatureSelectorScore
        #   or MRMFeatureSelectorQMIP directly in a loop:
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
        #       params = MRMFeatureSelector.SelectorParameters()
        #       params.nn_threshold = 4
        #       params.optimal_threshold = 0.5
        #
        #       # Run selection
        #       selector = MRMFeatureSelectorScore()
        #       selected = FeatureMap()
        #       selector.selectMRMFeature(features, selected, params)
        #
        #       # For iterative refinement, run again with stricter params
        #       params.optimal_threshold = 0.8
        #       final_selected = FeatureMap()
        #       selector.selectMRMFeature(selected, final_selected, params)

        # Note: Constructor is deleted in C++
        # Static methods are not wrapped due to nested type limitations
        pass
