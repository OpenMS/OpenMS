from Types cimport *
from Param cimport *
from AbsoluteQuantitationMethod cimport *
from AbsoluteQuantitationStandards cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from DefaultParamHandler cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitation(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #   Absolute quantitation using calibration curves and internal standards
        #
        #   This class supports absolute or relative quantitation for targeted workflows
        #   using Isotope Dilution Mass Spectrometry (IDMS). A transformation model is
        #   fitted using standards with known concentrations, then applied to quantify
        #   unknown samples.
        #
        #   Workflow:
        #   1. Set quantitation methods with setQuantMethods()
        #   2. Fit calibration curves with optimizeCalibrationCurves() or fitCalibration()
        #   3. Quantify unknowns with quantifyComponents()
        #
        #   Example usage:
        #
        #   .. code-block:: python
        #
        #       from pyopenms import *
        #
        #       # Set up quantitation method
        #       method = AbsoluteQuantitationMethod()
        #       method.setComponentName("glucose")
        #       method.setFeatureName("peak_apex_int")
        #       method.setISName("glucose_13C6")
        #       method.setTransformationModel("linear")
        #       method.setLLOQ(1.0)
        #       method.setULOQ(500.0)
        #
        #       aq = AbsoluteQuantitation()
        #       aq.setQuantMethods([method])
        #
        #       # Quantify unknowns
        #       unknowns = FeatureMap()
        #       FeatureXMLFile().load("unknowns.featureXML", unknowns)
        #       aq.quantifyComponents(unknowns)
        #
        #       # Results stored as metavalues
        #       for f in unknowns:
        #           if f.metaValueExists("calculated_concentration"):
        #               print(f.getMetaValue("PeptideRef"), f.getMetaValue("calculated_concentration"))

        AbsoluteQuantitation() except + nogil
        AbsoluteQuantitation(AbsoluteQuantitation &) except + nogil

        void setQuantMethods(libcpp_vector[AbsoluteQuantitationMethod] & quant_methods) except + nogil
            # wrap-doc:
            #   Set the quantitation methods for all components
            #
            #   :param quant_methods: Vector of AbsoluteQuantitationMethod objects

        libcpp_vector[AbsoluteQuantitationMethod] getQuantMethods() except + nogil
            # wrap-doc:
            #   Get the current quantitation methods
            #
            #   :returns: Vector of AbsoluteQuantitationMethod objects

        double calculateRatio(Feature & component_1, Feature & component_2, const String & feature_name) except + nogil
            # wrap-doc:
            #   Calculate the ratio between two features (e.g., analyte/IS)
            #
            #   :param component_1: First feature (numerator)
            #   :param component_2: Second feature (denominator, e.g., internal standard)
            #   :param feature_name: Name of the feature to use (e.g., "peak_apex_int")
            #   :returns: The ratio component_1 / component_2

        double calculateBias(double actual_concentration, double calculated_concentration) except + nogil
            # wrap-doc:
            #   Calculate the bias between actual and calculated concentration
            #
            #   Bias = (calculated - actual) / actual * 100%
            #
            #   :param actual_concentration: Known concentration
            #   :param calculated_concentration: Calculated concentration from calibration
            #   :returns: Bias as a percentage

        Param fitCalibration(
            libcpp_vector[AQS_featureConcentration] & component_concentrations,
            const String & feature_name,
            const String & transformation_model,
            Param transformation_model_params
        ) except + nogil
            # wrap-doc:
            #   Fit a calibration curve to standard concentrations
            #
            #   :param component_concentrations: Vector of featureConcentration with known values
            #   :param feature_name: Feature to use (e.g., "peak_apex_int")
            #   :param transformation_model: Model type (e.g., "linear", "b_spline")
            #   :param transformation_model_params: Initial model parameters
            #   :returns: Optimized model parameters

        void calculateBiasAndR(
            libcpp_vector[AQS_featureConcentration] & component_concentrations,
            const String & feature_name,
            const String & transformation_model,
            Param & transformation_model_params,
            libcpp_vector[double] & biases,
            double & correlation_coefficient
        ) except + nogil
            # wrap-doc:
            #   Calculate bias values and correlation coefficient for calibration
            #
            #   :param component_concentrations: Standards with known concentrations
            #   :param feature_name: Feature to use
            #   :param transformation_model: Model type
            #   :param transformation_model_params: Model parameters (updated in place)
            #   :param biases: Output vector of bias values for each standard
            #   :param correlation_coefficient: Output Pearson R correlation

        bool optimizeCalibrationCurveIterative(
            libcpp_vector[AQS_featureConcentration] & component_concentrations,
            const String & feature_name,
            const String & transformation_model,
            const Param & transformation_model_params,
            Param & optimized_params
        ) except + nogil
            # wrap-doc:
            #   Iteratively optimize calibration curve with outlier removal
            #
            #   Uses jackknife or residual-based outlier detection to iteratively
            #   remove points that degrade the calibration fit.
            #
            #   :param component_concentrations: Standards (outliers will be removed)
            #   :param feature_name: Feature to use
            #   :param transformation_model: Model type
            #   :param transformation_model_params: Initial model parameters
            #   :param optimized_params: Output optimized parameters
            #   :returns: True if optimization succeeded

        void optimizeCalibrationCurves(
            libcpp_map[String, libcpp_vector[AQS_featureConcentration]] & components_concentrations
        ) except + nogil
            # wrap-doc:
            #   Optimize calibration curves for all components
            #
            #   Fits and optimizes calibration curves for all components in the map.
            #   Results are stored in the internal quantitation methods.
            #
            #   :param components_concentrations: Map of component names to their standards

        void optimizeSingleCalibrationCurve(
            const String & component_name,
            libcpp_vector[AQS_featureConcentration] & component_concentrations
        ) except + nogil
            # wrap-doc:
            #   Optimize calibration curve for a single component
            #
            #   :param component_name: Name of the component
            #   :param component_concentrations: Standards for this component

        double applyCalibration(
            const Feature & component,
            const Feature & IS_component,
            const String & feature_name,
            const String & transformation_model,
            const Param & transformation_model_params
        ) except + nogil
            # wrap-doc:
            #   Apply calibration to calculate concentration
            #
            #   :param component: Analyte feature
            #   :param IS_component: Internal standard feature
            #   :param feature_name: Feature to use
            #   :param transformation_model: Model type
            #   :param transformation_model_params: Fitted model parameters
            #   :returns: Calculated concentration

        void quantifyComponents(FeatureMap & unknowns) except + nogil
            # wrap-doc:
            #   Quantify all components in a FeatureMap
            #
            #   Applies the calibration curves to calculate concentrations for all
            #   features. Results are stored as metavalues:
            #   - "calculated_concentration": The calculated concentration
            #   - "concentration_units": Units from the quantitation method
            #
            #   :param unknowns: FeatureMap to quantify (modified in place)
