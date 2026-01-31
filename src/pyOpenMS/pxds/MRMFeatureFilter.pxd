from Types cimport *
from Feature cimport *
from FeatureMap cimport *
from MRMFeatureQC cimport *
from DefaultParamHandler cimport *
from TargetedExperiment cimport *
from AbsoluteQuantitationMethod cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureFilter(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #   Flags or filters MRM features that do not pass QC criteria
        #
        #   This class provides comprehensive quality control filtering for MRM/SRM features.
        #   It can filter based on:
        #   - Retention time, intensity, and quality bounds
        #   - Ion ratios between transitions
        #   - %RSD from pooled QC samples
        #   - Background interference from blank samples
        #
        #   Example usage:
        #
        #   .. code-block:: python
        #
        #       from pyopenms import *
        #
        #       # Load features and transitions
        #       features = FeatureMap()
        #       FeatureXMLFile().load("features.featureXML", features)
        #       transitions = TargetedExperiment()
        #       TraMLFile().load("transitions.traML", transitions)
        #
        #       # Set up QC criteria
        #       qc = MRMFeatureQC()
        #       comp_qc = MRMFeatureQC.ComponentQCs()
        #       comp_qc.component_name = "my_transition"
        #       comp_qc.retention_time_l = 5.0
        #       comp_qc.retention_time_u = 7.0
        #       comp_qc.intensity_l = 1000.0
        #       qc.component_qcs.append(comp_qc)
        #
        #       # Apply filtering
        #       filter = MRMFeatureFilter()
        #       filter.FilterFeatureMap(features, qc, transitions)

        MRMFeatureFilter() except + nogil
        MRMFeatureFilter(MRMFeatureFilter &) except + nogil

        void FilterFeatureMap(FeatureMap & features, MRMFeatureQC & filter_criteria, TargetedExperiment & transitions) except + nogil
            # wrap-doc:
            #   Flags or filters features and subordinates in a FeatureMap
            #
            #   :param features: FeatureMap to flag or filter
            #   :param filter_criteria: MRMFeatureQC class defining QC parameters
            #   :param transitions: Transitions from a TargetedExperiment

        void FilterFeatureMapPercRSD(FeatureMap & features, MRMFeatureQC & filter_criteria, MRMFeatureQC & filter_values) except + nogil
            # wrap-doc:
            #   Filter features based on percent relative standard deviation (%RSD)
            #
            #   Uses QC samples to calculate %RSD and filters features exceeding thresholds.
            #
            #   :param features: FeatureMap to filter
            #   :param filter_criteria: MRMFeatureQC with %RSD thresholds
            #   :param filter_values: MRMFeatureQC with calculated %RSD values

        void FilterFeatureMapBackgroundInterference(FeatureMap & features, MRMFeatureQC & filter_criteria, MRMFeatureQC & filter_values) except + nogil
            # wrap-doc:
            #   Filter features based on background interference from blank samples
            #
            #   :param features: FeatureMap to filter
            #   :param filter_criteria: MRMFeatureQC with interference thresholds
            #   :param filter_values: MRMFeatureQC with calculated interference values

        void EstimateDefaultMRMFeatureQCValues(
            libcpp_vector[FeatureMap] & samples,
            MRMFeatureQC & filter_template,
            TargetedExperiment & transitions,
            bool init_template_values
        ) except + nogil
            # wrap-doc:
            #   Estimate default QC values from a set of standard samples
            #
            #   Analyzes standard samples to automatically determine appropriate QC bounds.
            #
            #   :param samples: Vector of FeatureMaps from standard samples
            #   :param filter_template: MRMFeatureQC to populate with estimated values
            #   :param transitions: TargetedExperiment with transition information
            #   :param init_template_values: If true, initialize template values

        void TransferLLOQAndULOQToCalculatedConcentrationBounds(
            libcpp_vector[AbsoluteQuantitationMethod] & quantitation_method,
            MRMFeatureQC & filter_template
        ) except + nogil
            # wrap-doc:
            #   Transfer LLOQ/ULOQ from quantitation methods to QC bounds
            #
            #   Copies the lower and upper limits of quantitation from AbsoluteQuantitationMethod
            #   objects to the corresponding concentration bounds in the MRMFeatureQC template.
            #
            #   :param quantitation_method: Vector of AbsoluteQuantitationMethod with LOQ values
            #   :param filter_template: MRMFeatureQC to update with concentration bounds

        void EstimatePercRSD(
            libcpp_vector[FeatureMap] & samples,
            MRMFeatureQC & filter_template,
            TargetedExperiment & transitions
        ) except + nogil
            # wrap-doc:
            #   Estimate percent relative standard deviation from QC samples
            #
            #   Calculates %RSD across multiple QC sample injections.
            #
            #   :param samples: Vector of FeatureMaps from pooled QC samples
            #   :param filter_template: MRMFeatureQC to populate with %RSD values
            #   :param transitions: TargetedExperiment with transition information

        void EstimateBackgroundInterferences(
            libcpp_vector[FeatureMap] & samples,
            MRMFeatureQC & filter_template,
            TargetedExperiment & transitions
        ) except + nogil
            # wrap-doc:
            #   Estimate background interference from blank samples
            #
            #   Analyzes blank samples to determine background signal levels.
            #
            #   :param samples: Vector of FeatureMaps from blank samples
            #   :param filter_template: MRMFeatureQC to populate with interference values
            #   :param transitions: TargetedExperiment with transition information

        double calculateIonRatio(const Feature & component_1, const Feature & component_2, const String & feature_name) except + nogil
            # wrap-doc:
            #   Calculate the ion ratio between two features
            #
            #   :param component_1: First feature (e.g., quantifier transition)
            #   :param component_2: Second feature (e.g., qualifier transition)
            #   :param feature_name: Name of the feature to use for ratio (e.g., "peak_apex_int")
            #   :returns: The ratio of component_1 / component_2

        double calculateRTDifference(Feature & component_1, Feature & component_2) except + nogil
            # wrap-doc:
            #   Calculate the retention time difference between two features
            #
            #   :param component_1: First feature
            #   :param component_2: Second feature
            #   :returns: The absolute RT difference in seconds

        double calculateResolution(Feature & component_1, Feature & component_2) except + nogil
            # wrap-doc:
            #   Calculate the chromatographic resolution between two features
            #
            #   Resolution is calculated as: 2 * |RT2 - RT1| / (W1 + W2)
            #   where W is the peak width at base (estimated from FWHM if available).
            #
            #   :param component_1: First feature
            #   :param component_2: Second feature
            #   :returns: The resolution value (0 if peak widths are not available)
