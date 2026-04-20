// pyOpenMS nanobind bindings
// Domain: misc

#include "all_casters.h"
#include "nanobind_ms_data_consumer.h"
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h>
#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmIdentity.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>
#include <OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/ANALYSIS/XLMS/OpenPepXLAlgorithm.h>
#include <OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
#include <OpenMS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
#include <OpenMS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/FEATUREFINDER/LevMarqFitter1D.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataAggregatingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataStoringConsumer.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/GNPSMGFFile.h>
#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MRMFeaturePickerFile.h>
#include <OpenMS/FORMAT/MRMFeatureQCFile.h>
#include <OpenMS/FORMAT/MS2File.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/PEFFFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFileMascot.h>
#include <OpenMS/FORMAT/ProtXMLFile.h>
#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MISC/EmgGradientDescent.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ML/CLUSTERING/AverageLinkage.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/PROCESSING/SCALING/RankScaler.h>
#include <OpenMS/PROCESSING/SCALING/SqrtScaler.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>
#include <OpenMS/PROCESSING/SMOOTHING/LowessSmoothing.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/SYSTEM/BuildInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>
#include <iomanip>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/set.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>
#include <sstream>
#include "binding_utils.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_pyopenms_misc, m) {
    m.doc() = "pyOpenMS misc bindings";

    // -----------------------------------------------------------------------
    // BinaryDataArray
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Interfaces::BinaryDataArray>(m, "BinaryDataArray", "OpenMS class BinaryDataArray")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Interfaces::BinaryDataArray &>())
        .def("__copy__", [](const OpenMS::Interfaces::BinaryDataArray& self) { return OpenMS::Interfaces::BinaryDataArray(self); })
        .def("__deepcopy__", [](const OpenMS::Interfaces::BinaryDataArray& self, nb::dict) { return OpenMS::Interfaces::BinaryDataArray(self); }, "memo"_a)
        .def_rw("data", &OpenMS::Interfaces::BinaryDataArray::data)
        ;

    // -----------------------------------------------------------------------
    // Chromatogram
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Interfaces::Chromatogram>(m, "Chromatogram", "OpenMS class Chromatogram")
        .def(nb::init<>())
        .def_rw("defaultArrayLength", &OpenMS::Interfaces::Chromatogram::defaultArrayLength)

        .def("setTimeArray", [](OpenMS::Interfaces::Chromatogram& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setTimeArray(arr);
        }, "data"_a, "Set time array from list")

        .def("setIntensityArray", [](OpenMS::Interfaces::Chromatogram& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")

        .def("getTimeArray", [](const OpenMS::Interfaces::Chromatogram& self) {
            auto arr = self.getTimeArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get time array")

        .def("getIntensityArray", [](const OpenMS::Interfaces::Chromatogram& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array")
        ;

    // -----------------------------------------------------------------------
    // DefaultParamHandler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DefaultParamHandler>(m, "DefaultParamHandler", "A base class for all classes handling default parameters")
        .def(nb::init<const OpenMS::DefaultParamHandler &>())
        .def(nb::init<OpenMS::String>())
        .def("__copy__", [](const OpenMS::DefaultParamHandler& self) { return OpenMS::DefaultParamHandler(self); })
        .def("__deepcopy__", [](const OpenMS::DefaultParamHandler& self, nb::dict) { return OpenMS::DefaultParamHandler(self); }, "memo"_a)
        .def("setParameters", [](OpenMS::DefaultParamHandler& self, const OpenMS::Param& param) { return self.setParameters(param); }, "param"_a, "Sets the parameters")
        .def("getParameters", [](const OpenMS::DefaultParamHandler& self) -> const OpenMS::Param & { return self.getParameters(); }, nb::rv_policy::reference_internal, "Returns the parameters")
        .def("getDefaults", [](const OpenMS::DefaultParamHandler& self) -> const OpenMS::Param & { return self.getDefaults(); }, nb::rv_policy::reference_internal, "Returns the default parameters")
        .def("getName", [](const OpenMS::DefaultParamHandler& self) { return self.getName(); }, "Returns the name")
        .def("setName", [](OpenMS::DefaultParamHandler& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name")
        .def("getSubsections", [](const OpenMS::DefaultParamHandler& self) -> const std::vector<OpenMS::String> & { return self.getSubsections(); }, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // AScore
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AScore, OpenMS::DefaultParamHandler>(m, "AScore", 
        R"doc(
Implementation of the Ascore For a given peptide sequence and its
MS/MS spectrum it identifies the most probable
phosphorylation-site(s). For each phosphorylation site a probability
score is calculated. The algorithm is implemented according to
Beausoleil et al. (Nat. Biotechnol. 2006)
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("compute", [](OpenMS::AScore& self, const OpenMS::PeptideHit& hit, OpenMS::MSSpectrum& real_spectrum) { return self.compute(hit, real_spectrum); }, "hit"_a, "real_spectrum"_a)
        ;

    // -----------------------------------------------------------------------
    // AbsoluteQuantitation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitation, OpenMS::DefaultParamHandler>(m, "AbsoluteQuantitation", 
        R"doc(
DefaultParamHandler

Absolute quantitation using calibration curves and internal standards
This class supports absolute or relative quantitation for targeted workflows
using Isotope Dilution Mass Spectrometry (IDMS). A transformation model is
fitted using standards with known concentrations, then applied to quantify
unknown samples.
Workflow:
1. Set quantitation methods with setQuantMethods()
2. Fit calibration curves with optimizeCalibrationCurves() or fitCalibration()
3. Quantify unknowns with quantifyComponents()
Example usage:
.. code-block:: python
from pyopenms import *
# Set up quantitation method
method = AbsoluteQuantitationMethod()
method.setComponentName("glucose")
method.setFeatureName("peak_apex_int")
method.setISName("glucose_13C6")
method.setTransformationModel("linear")
method.setLLOQ(1.0)
method.setULOQ(500.0)
aq = AbsoluteQuantitation()
aq.setQuantMethods([method])
# Quantify unknowns
unknowns = FeatureMap()
FeatureXMLFile().load("unknowns.featureXML", unknowns)
aq.quantifyComponents(unknowns)
# Results stored as metavalues
for f in unknowns:
if f.metaValueExists("calculated_concentration"):
print(f.getMetaValue("PeptideRef"), f.getMetaValue("calculated_concentration"))
)doc")
        .def(nb::init<>())
        .def("calculateRatio", [](OpenMS::AbsoluteQuantitation& self, const OpenMS::Feature& component_1, const OpenMS::Feature& component_2, const OpenMS::String& feature_name) { return self.calculateRatio(component_1, component_2, feature_name); }, "component_1"_a, "component_2"_a, "feature_name"_a, 
            R"doc(
Get the current quantitation methods
:returns: Vector of AbsoluteQuantitationMethod objects
)doc")
        .def("calculateBias", [](OpenMS::AbsoluteQuantitation& self, const double& actual_concentration, const double& calculated_concentration) { return self.calculateBias(actual_concentration, calculated_concentration); }, "actual_concentration"_a, "calculated_concentration"_a, 
            R"doc(
Calculate the ratio between two features (e.g., analyte/IS)
:param component_1: First feature (numerator)
:param component_2: Second feature (denominator, e.g., internal standard)
:param feature_name: Name of the feature to use (e.g., "peak_apex_int")
:returns: The ratio component_1 / component_2
)doc")
        .def("optimizeCalibrationCurveIterative", [](OpenMS::AbsoluteQuantitation& self, std::vector<OpenMS::AbsoluteQuantitationStandards::featureConcentration> component_concentrations, const OpenMS::String& feature_name, const OpenMS::String& transformation_model, const OpenMS::Param& transformation_model_params, OpenMS::Param& optimized_params) { bool result = self.optimizeCalibrationCurveIterative(component_concentrations, feature_name, transformation_model, transformation_model_params, optimized_params); return nb::make_tuple(result, component_concentrations); }, "component_concentrations"_a, "feature_name"_a, "transformation_model"_a, "transformation_model_params"_a, "optimized_params"_a,
            R"doc(
Calculate bias values and correlation coefficient for calibration
:param component_concentrations: Standards with known concentrations
:param feature_name: Feature to use
:param transformation_model: Model type
:param transformation_model_params: Model parameters
:param optimized_params: Output optimized parameters (modified in-place)
:returns: Tuple of (success_bool, updated_component_concentrations)
)doc")
        .def("optimizeSingleCalibrationCurve", [](OpenMS::AbsoluteQuantitation& self, const OpenMS::String& component_name, std::vector<OpenMS::AbsoluteQuantitationStandards::featureConcentration> component_concentrations) { self.optimizeSingleCalibrationCurve(component_name, component_concentrations); return component_concentrations; }, "component_name"_a, "component_concentrations"_a,
            R"doc(
Iteratively optimize calibration curve with outlier removal
Uses jackknife or residual-based outlier detection to iteratively
remove points that degrade the calibration fit.
:param component_concentrations: Standards (outliers will be removed)
:returns: Updated component_concentrations after optimization
Note: optimizeCalibrationCurves taking map<String, vector<featureConcentration>>
is not wrapped due to complex nested container type. Use optimizeSingleCalibrationCurve instead.
)doc")
        .def("applyCalibration", [](OpenMS::AbsoluteQuantitation& self, const OpenMS::Feature& component, const OpenMS::Feature& IS_component, const OpenMS::String& feature_name, const OpenMS::String& transformation_model, const OpenMS::Param& transformation_model_params) { return self.applyCalibration(component, IS_component, feature_name, transformation_model, transformation_model_params); }, "component"_a, "IS_component"_a, "feature_name"_a, "transformation_model"_a, "transformation_model_params"_a, 
            R"doc(
Optimize calibration curve for a single component
:param component_name: Name of the component
:param component_concentrations: Standards for this component
)doc")
        .def("quantifyComponents", [](OpenMS::AbsoluteQuantitation& self, OpenMS::FeatureMap& unknowns) { return self.quantifyComponents(unknowns); }, "unknowns"_a, 
            R"doc(
Apply calibration to calculate concentration
:param component: Analyte feature
:param IS_component: Internal standard feature
:param feature_name: Feature to use
:param transformation_model: Model type
:param transformation_model_params: Fitted model parameters
:returns: Calculated concentration
)doc")

        .def("setQuantMethods", [](OpenMS::AbsoluteQuantitation& self, const std::vector<OpenMS::AbsoluteQuantitationMethod>& quant_methods) {
            std::vector<OpenMS::AbsoluteQuantitationMethod> methods_copy(quant_methods);
            self.setQuantMethods(methods_copy);
        }, "quant_methods"_a, "Set the quantitation methods")

        .def("getQuantMethods", [](OpenMS::AbsoluteQuantitation& self) {
            return self.getQuantMethods();
        }, "Get the quantitation methods")
        .def("fitCalibration", [](OpenMS::AbsoluteQuantitation& self, const std::vector<OpenMS::AbsoluteQuantitationStandards::featureConcentration>& component_concentrations, const OpenMS::String& feature_name, const OpenMS::String& transformation_model, const OpenMS::Param& transformation_model_params) { return self.fitCalibration(component_concentrations, feature_name, transformation_model, transformation_model_params); }, "component_concentrations"_a, "feature_name"_a, "transformation_model"_a, "transformation_model_params"_a, "Fit a calibration curve")
        .def("calculateBiasAndR", [](OpenMS::AbsoluteQuantitation& self, const std::vector<OpenMS::AbsoluteQuantitationStandards::featureConcentration>& component_concentrations, const OpenMS::String& feature_name, const OpenMS::String& transformation_model, const OpenMS::Param& transformation_model_params) {
            std::vector<double> biases;
            double correlation_coefficient;
            self.calculateBiasAndR(component_concentrations, feature_name, transformation_model, transformation_model_params, biases, correlation_coefficient);
            return nb::make_tuple(biases, correlation_coefficient);
        }, "component_concentrations"_a, "feature_name"_a, "transformation_model"_a, "transformation_model_params"_a, "Calculate bias values and R for calibration")
        ;

    // -----------------------------------------------------------------------
    // Biosaur2Algorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Biosaur2Algorithm, OpenMS::DefaultParamHandler>(m, "Biosaur2Algorithm", 
        R"doc(
DefaultParamHandler

C++ implementation of the Biosaur2 feature detection workflow.
)doc")
        .def(nb::init<>())
        .def("setMSData", [](OpenMS::Biosaur2Algorithm& self, const OpenMS::MSExperiment& ms_data) { return self.setMSData(ms_data); }, "ms_data"_a, "Set the MS data used for feature detection (copy version)")
        .def("setMSData", [](OpenMS::Biosaur2Algorithm& self, OpenMS::MSExperiment& ms_data) { return self.setMSData(ms_data); }, "ms_data"_a, "Set the MS data used for feature detection (copy version)")
        .def("getMSData", [](OpenMS::Biosaur2Algorithm& self) -> OpenMS::MSExperiment & { return self.getMSData(); }, nb::rv_policy::reference_internal, "Get non-const reference to MS data")
        .def("run", [](OpenMS::Biosaur2Algorithm& self, OpenMS::FeatureMap& feature_map) { nb::gil_scoped_release release; return self.run(feature_map); }, "feature_map"_a, "Run the algorithm storing only the resulting features")
        .def("run", [](OpenMS::Biosaur2Algorithm& self, OpenMS::FeatureMap& feature_map, bool /*return_details*/) { std::vector<OpenMS::Biosaur2Algorithm::Hill> hills; std::vector<OpenMS::Biosaur2Algorithm::PeptideFeature> peptide_features; { nb::gil_scoped_release release; self.run(feature_map, hills, peptide_features); } return nb::make_tuple(hills, peptide_features); }, "feature_map"_a, "return_details"_a, "Run the algorithm returning (hills, peptide_features) in addition to modifying feature_map in-place")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithm, OpenMS::DefaultParamHandler>(m, "ConsensusIDAlgorithm", 
        R"doc(
Abstract base class for all ConsensusID algorithms (that calculate a
consensus from multiple ID runs)
DefaultParamHandler
)doc")
        .def("apply", [](OpenMS::ConsensusIDAlgorithm& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a = 0, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        .def("apply", [](OpenMS::ConsensusIDAlgorithm& self, OpenMS::PeptideIdentificationList& ids, size_t number_of_runs) { return self.apply(ids, number_of_runs); }, "ids"_a, "number_of_runs"_a = 0, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmIdentity
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmIdentity, OpenMS::ConsensusIDAlgorithm>(m, "ConsensusIDAlgorithmIdentity", 
        R"doc(
Abstract base class for ConsensusID algorithms that compare only
identical sequences
ConsensusIDAlgorithm
)doc")
        .def("apply", [](OpenMS::ConsensusIDAlgorithmIdentity& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmAverage
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmAverage, OpenMS::ConsensusIDAlgorithmIdentity>(m, "ConsensusIDAlgorithmAverage", 
        R"doc(
Calculates a consensus from multiple ID runs by averaging the search
scores
ConsensusIDAlgorithmIdentity
)doc")
        .def(nb::init<>())
        .def("apply", [](OpenMS::ConsensusIDAlgorithmAverage& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmBest
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmBest, OpenMS::ConsensusIDAlgorithmIdentity>(m, "ConsensusIDAlgorithmBest", 
        R"doc(
Calculates a consensus from multiple ID runs by taking the best
search score
ConsensusIDAlgorithmIdentity
)doc")
        .def(nb::init<>())
        .def("apply", [](OpenMS::ConsensusIDAlgorithmBest& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmRanks
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmRanks, OpenMS::ConsensusIDAlgorithmIdentity>(m, "ConsensusIDAlgorithmRanks", 
        R"doc(
Calculates a consensus from multiple ID runs based on the ranks of
the search hits
ConsensusIDAlgorithmIdentity
)doc")
        .def(nb::init<>())
        .def("apply", [](OpenMS::ConsensusIDAlgorithmRanks& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmSimilarity
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmSimilarity, OpenMS::ConsensusIDAlgorithm>(m, "ConsensusIDAlgorithmSimilarity", 
        R"doc(
Abstract base class for ConsensusID algorithms that take peptide
similarity into account
ConsensusIDAlgorithm
)doc")
        .def("apply", [](OpenMS::ConsensusIDAlgorithmSimilarity& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmPEPIons
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmPEPIons, OpenMS::ConsensusIDAlgorithmSimilarity>(m, "ConsensusIDAlgorithmPEPIons", 
        R"doc(
Calculates a consensus from multiple ID runs based on PEPs and shared
ions
ConsensusIDAlgorithmSimilarity
)doc")
        .def(nb::init<>())
        .def("apply", [](OpenMS::ConsensusIDAlgorithmPEPIons& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmPEPMatrix
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmPEPMatrix, OpenMS::ConsensusIDAlgorithmSimilarity>(m, "ConsensusIDAlgorithmPEPMatrix", 
        R"doc(
Calculates a consensus from multiple ID runs based on PEPs and
sequence similarities
ConsensusIDAlgorithmSimilarity
)doc")
        .def(nb::init<>())
        .def("apply", [](OpenMS::ConsensusIDAlgorithmPEPMatrix& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ConsensusIDAlgorithmWorst
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusIDAlgorithmWorst, OpenMS::ConsensusIDAlgorithmIdentity>(m, "ConsensusIDAlgorithmWorst", 
        R"doc(
Calculates a consensus from multiple ID runs by taking the worst
search score (conservative approach)
ConsensusIDAlgorithmIdentity
)doc")
        .def(nb::init<>())
        .def("apply", [](OpenMS::ConsensusIDAlgorithmWorst& self, OpenMS::PeptideIdentificationList& ids, const std::map<OpenMS::String, OpenMS::String>& se_info, size_t number_of_runs) { return self.apply(ids, se_info, number_of_runs); }, "ids"_a, "se_info"_a, "number_of_runs"_a, "Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature")
        ;

    // -----------------------------------------------------------------------
    // ElutionModelFitter
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ElutionModelFitter, OpenMS::DefaultParamHandler>(m, "ElutionModelFitter", 
        R"doc(
Helper class for fitting elution models to features
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("fitElutionModels", [](OpenMS::ElutionModelFitter& self, OpenMS::FeatureMap& features) { return self.fitElutionModels(features); }, "features"_a, "Fit models of elution profiles to all features (and validate them)")
        ;

    // -----------------------------------------------------------------------
    // EmgGradientDescent
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::EmgGradientDescent, OpenMS::DefaultParamHandler>(m, "EmgGradientDescent", 
        R"doc(
Fit peaks to an Exponentially Modified Gaussian (EMG) model using
gradient descent
DefaultParamHandler
)doc")
        .def(nb::init<>())

        .def("getDefaultParameters", [](OpenMS::EmgGradientDescent& self, OpenMS::Param& params) {
            self.getDefaultParameters(params);
        }, "params"_a, "Get default parameters (fills output param)")
        .def("getDefaultParameters", [](OpenMS::EmgGradientDescent& self) {
            OpenMS::Param params;
            self.getDefaultParameters(params);
            return params;
        }, "Get default parameters (returns Param)")
        .def("fitEMGPeakModel", [](const OpenMS::EmgGradientDescent& self, const OpenMS::MSChromatogram& input_peak, OpenMS::MSChromatogram& output_peak, double left_pos, double right_pos) { self.fitEMGPeakModel(input_peak, output_peak, left_pos, right_pos); }, "input_peak"_a, "output_peak"_a, "left_pos"_a = 0.0, "right_pos"_a = 0.0, "Fit EMG peak model to chromatographic peak")
        ;

    // -----------------------------------------------------------------------
    // FIAMSDataProcessor
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FIAMSDataProcessor, OpenMS::DefaultParamHandler>(m, "FIAMSDataProcessor", 
        R"doc(
Processes FIA-MS data for metabolite identification
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FIAMSDataProcessor &>())
        .def("__copy__", [](const OpenMS::FIAMSDataProcessor& self) { return OpenMS::FIAMSDataProcessor(self); })
        .def("__deepcopy__", [](const OpenMS::FIAMSDataProcessor& self, nb::dict) { return OpenMS::FIAMSDataProcessor(self); }, "memo"_a)
        .def("run", [](OpenMS::FIAMSDataProcessor& self, const OpenMS::MSExperiment& experiment, float n_seconds, OpenMS::MzTab& output, bool load_cached_spectrum) { return self.run(experiment, n_seconds, output, load_cached_spectrum); }, "experiment"_a, "n_seconds"_a, "output"_a, "load_cached_spectrum"_a = true)
        .def("extractPeaks", [](OpenMS::FIAMSDataProcessor& self, const OpenMS::MSSpectrum& input) { return self.extractPeaks(input); }, "input"_a, 
            R"doc(
Run the full analysis for the experiment for the given time interval\n
The workflow steps are:
- the time axis of the experiment is cut to the interval from 0 to n_seconds
- the spectra are summed into one along the time axis with the bin size determined by mz and instrument resolution
- data is smoothed by applying the Savitzky-Golay filter
- peaks are picked
- the accurate mass search for all the picked peaks is performed
The intermediate summed spectra and picked peaks can be saved to the filesystem.
Also, the results of the accurate mass search and the signal-to-noise information
of the resulting spectrum is saved.
:param experiment: Input MSExperiment
:param n_seconds: Input number of seconds
:param load_cached_spectrum: Load the cached picked spectrum if exists
:param output: Output of the accurate mass search results
:return: A boolean indicating if the picked spectrum was loaded from the cached file
)doc")
        .def("convertToFeatureMap", [](OpenMS::FIAMSDataProcessor& self, const OpenMS::MSSpectrum& input) { return self.convertToFeatureMap(input); }, "input"_a, 
            R"doc(
Pick peaks from the summed spectrum
:param input: Input vector of spectra
:return: A spectrum with picked peaks
)doc")
        .def("trackNoise", [](OpenMS::FIAMSDataProcessor& self, const OpenMS::MSSpectrum& input) { return self.trackNoise(input); }, "input"_a, 
            R"doc(
Convert a spectrum to a feature map with the corresponding polarity\n
Applies `SavitzkyGolayFilter` and `PeakPickerHiRes`
:param input: Input a picked spectrum
:return: A feature map with the peaks converted to features and polarity from the parameters
)doc")
        ;

    // -----------------------------------------------------------------------
    // FalseDiscoveryRate
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FalseDiscoveryRate, OpenMS::DefaultParamHandler>(m, "FalseDiscoveryRate", 
        R"doc(
Calculates false discovery rates (FDR) from identifications
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("apply", [](const OpenMS::FalseDiscoveryRate& self, OpenMS::PeptideIdentificationList& fwd_ids, OpenMS::PeptideIdentificationList& rev_ids) { return self.apply(fwd_ids, rev_ids); }, "fwd_ids"_a, "rev_ids"_a)
        .def("apply", [](const OpenMS::FalseDiscoveryRate& self, OpenMS::PeptideIdentificationList& id, bool annotate_peptide_fdr) { return self.apply(id, annotate_peptide_fdr); }, "id"_a, "annotate_peptide_fdr"_a = false)
        .def("apply", [](const OpenMS::FalseDiscoveryRate& self, std::vector<OpenMS::ProteinIdentification> fwd_ids, std::vector<OpenMS::ProteinIdentification> rev_ids) {
            self.apply(fwd_ids, rev_ids);
            return nb::make_tuple(fwd_ids, rev_ids);
        }, "fwd_ids"_a, "rev_ids"_a)
        .def("apply", [](const OpenMS::FalseDiscoveryRate& self, std::vector<OpenMS::ProteinIdentification> ids) {
            self.apply(ids);
            return ids;
        }, "ids"_a)
        .def("applyEvaluateProteinIDs", [](const OpenMS::FalseDiscoveryRate& self, const std::vector<OpenMS::ProteinIdentification>& ids, double pepCutoff, unsigned int fpCutoff, double diffWeight) { return self.applyEvaluateProteinIDs(ids, pepCutoff, fpCutoff, diffWeight); }, "ids"_a, "pepCutoff"_a = 1.0, "fpCutoff"_a = 50, "diffWeight"_a = 0.2)
        .def("applyEvaluateProteinIDs", [](const OpenMS::FalseDiscoveryRate& self, const OpenMS::ProteinIdentification& ids, double pepCutoff, unsigned int fpCutoff, double diffWeight) { return self.applyEvaluateProteinIDs(ids, pepCutoff, fpCutoff, diffWeight); }, "ids"_a, "pepCutoff"_a = 1.0, "fpCutoff"_a = 50, "diffWeight"_a = 0.2)
        .def("applyPickedProteinFDR", [](OpenMS::FalseDiscoveryRate& self, OpenMS::ProteinIdentification& id, OpenMS::String decoy_string, bool prefix, bool groups_too) { return self.applyPickedProteinFDR(id, decoy_string, prefix, groups_too); }, "id"_a, "decoy_string"_a = "", "prefix"_a = true, "groups_too"_a = true)
        .def("applyBasic", [](OpenMS::FalseDiscoveryRate& self, OpenMS::PeptideIdentificationList& ids, bool higher_score_better, int charge, OpenMS::String identifier, bool only_best_per_pep) { self.applyBasic(ids, higher_score_better, charge, identifier, only_best_per_pep); }, "ids"_a, "higher_score_better"_a, "charge"_a = 0, "identifier"_a = "", "only_best_per_pep"_a = false, "Applies basic FDR calculation")
        .def("applyBasic", [](OpenMS::FalseDiscoveryRate& self, const std::vector<OpenMS::ProteinIdentification>& run_info, OpenMS::PeptideIdentificationList& ids) { self.applyBasic(run_info, ids); }, "run_info"_a, "ids"_a, "Applies basic FDR calculation using run info")
        .def("applyBasic", [](OpenMS::FalseDiscoveryRate& self, OpenMS::ConsensusMap& cmap, bool use_unassigned_peptides) { self.applyBasic(cmap, use_unassigned_peptides); }, "cmap"_a, "use_unassigned_peptides"_a = true, "Applies basic FDR calculation on ConsensusMap")
        .def("applyBasic", [](OpenMS::FalseDiscoveryRate& self, OpenMS::ProteinIdentification& id, bool groups_too) { self.applyBasic(id, groups_too); }, "id"_a, "groups_too"_a = true, "Applies basic FDR calculation on protein level")
        .def("applyBasicPeptideLevel", [](OpenMS::FalseDiscoveryRate& self, OpenMS::PeptideIdentificationList& ids) { self.applyBasicPeptideLevel(ids); }, "ids"_a, "Applies basic peptide-level FDR calculation")
        .def("applyBasicPeptideLevel", [](OpenMS::FalseDiscoveryRate& self, OpenMS::ConsensusMap& cmap, bool use_unassigned_peptides) { self.applyBasicPeptideLevel(cmap, use_unassigned_peptides); }, "cmap"_a, "use_unassigned_peptides"_a = true, "Applies basic peptide-level FDR calculation on ConsensusMap")
        .def("applyEstimated", [](const OpenMS::FalseDiscoveryRate& self, std::vector<OpenMS::ProteinIdentification> ids) { self.applyEstimated(ids); return ids; }, "ids"_a, "Applies estimated FDR calculation on protein IDs")
        .def("rocN", [](const OpenMS::FalseDiscoveryRate& self, const OpenMS::PeptideIdentificationList& ids, size_t fp_cutoff) { return self.rocN(ids, fp_cutoff); }, "ids"_a, "fp_cutoff"_a, "Calculates the ROC-N value (AUC) for peptide IDs")
        .def("rocN", [](const OpenMS::FalseDiscoveryRate& self, const OpenMS::PeptideIdentificationList& ids, size_t fp_cutoff, const OpenMS::String& identifier) { return self.rocN(ids, fp_cutoff, identifier); }, "ids"_a, "fp_cutoff"_a, "identifier"_a, "Calculates the ROC-N value for a specific identifier")
        ;

    // -----------------------------------------------------------------------
    // FeatureDistance
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureDistance, OpenMS::DefaultParamHandler>(m, "FeatureDistance", 
        R"doc(
A functor class for the calculation of distances between features or
consensus features
DefaultParamHandler
)doc")
        .def(nb::init<double, bool>())
        ;

    // -----------------------------------------------------------------------
    // FeatureFinderAlgorithmMetaboIdent
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureFinderAlgorithmMetaboIdent, OpenMS::DefaultParamHandler>(m, "FeatureFinderAlgorithmMetaboIdent", 
        R"doc(
DefaultParamHandler

Perform targeted feature extraction of compounds provided as table and stores them in features
The algorithms detects quantitative features in MS1 data for a list of targets, typically small molecule/metabolite identifications
Internally, it uses algorithms for targeted data analysis from the OpenSWATH pipeline
In the simplest case, only CompoundName, SumFormula, Charge and RetentionTime need to be given, all other values may be zero
Every combination of compound (mass), RT and charge defines one target for feature detection
Output:
The main output is a feature map of detected features, with annotations in meta data entries
Additional outputs are the extracted chromatograms/peak groups, the assay in TraML compatible format, and transformations
that contain the error between provided and observed peaks
Usage:
.. code-block:: python
exp = MSExperiment()
MzMLFile().load(path_to_file, exp)
ff = FeatureFinderAlgorithmMetaboIdent()
ff.setMSData(exp)
fm = FeatureMap() # detected features will be stored here
library = []
# fill library with compounds: FeatureFinderMetaboIdentCompound(name, formula, mass, [charges] [RTs_in_sec], [RT_ranges], [isotope distributions])
# e.g. FeatureFinderMetaboIdentCompound('glucose','C6H12O6', 0.0, [-1], [123.4], [0.0], [0.0])
params = ff.getParameters() # optional!
params[param_name] = new_value # e.g. params[b'extract:n_isotopes'] = 3
ff.setParameters(params)
ff.run(library, fm, path_to_file)
)doc")
        .def(nb::init<>())
        .def("getMSData", [](OpenMS::FeatureFinderAlgorithmMetaboIdent& self) -> OpenMS::MSExperiment & { return self.getMSData(); }, nb::rv_policy::reference_internal, "Returns spectra")
        .def("setMSData", [](OpenMS::FeatureFinderAlgorithmMetaboIdent& self, const OpenMS::MSExperiment& m) { return self.setMSData(m); }, "m"_a, "Sets spectra")
        .def("setMSData", [](OpenMS::FeatureFinderAlgorithmMetaboIdent& self, OpenMS::MSExperiment& m) { return self.setMSData(m); }, "m"_a, "Sets spectra")
        .def("getChromatograms", [](OpenMS::FeatureFinderAlgorithmMetaboIdent& self) -> OpenMS::MSExperiment & { return self.getChromatograms(); }, nb::rv_policy::reference_internal, "Retrieves chromatograms (empty if run was not executed)")
        .def("getLibrary", [](const OpenMS::FeatureFinderAlgorithmMetaboIdent& self) -> const OpenMS::TargetedExperiment & { return self.getLibrary(); }, nb::rv_policy::reference_internal, "Retrieves the assay library (e.g., to store as TraML, empty if run was not executed)")
        .def("getTransformations", [](const OpenMS::FeatureFinderAlgorithmMetaboIdent& self) -> const OpenMS::TransformationDescription & { return self.getTransformations(); }, nb::rv_policy::reference_internal, "Retrieves deviations between provided coordinates and extacted ones (e.g., to store as TrafoXML or for plotting)")
        .def("getNShared", [](const OpenMS::FeatureFinderAlgorithmMetaboIdent& self) { return self.getNShared(); }, "Retrieves number of features with shared identifications")
        .def("run", [](OpenMS::FeatureFinderAlgorithmMetaboIdent& self,
                       const std::vector<OpenMS::FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound>& metaboIdentTable,
                       OpenMS::FeatureMap& features,
                       const OpenMS::String& spectra_file) {
            self.run(metaboIdentTable, features, spectra_file);
        }, "metaboIdentTable"_a, "features"_a, "spectra_file"_a = "",
            R"doc(Perform targeted feature extraction of compounds and store them in features.

If spectra_file is provided it will be used as a fall-back to setPrimaryMSRunPath
in the feature map in case a proper primaryMSRunPath is not annotated in the MSExperiment.
If there are no MS1 scans in the MSData, features will be returned unchanged.

:param metaboIdentTable: List of FeatureFinderMetaboIdentCompound objects defining targets
:param features: FeatureMap to store detected features (modified in-place)
:param spectra_file: Optional path to spectra file for annotation
)doc")
        ;

    // -----------------------------------------------------------------------
    // FeatureFinderIdentificationAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureFinderIdentificationAlgorithm, OpenMS::DefaultParamHandler>(m, "FeatureFinderIdentificationAlgorithm", 
        R"doc(
DefaultParamHandler

Algorithm class for FeatureFinderIdentification
External IDs (peptides_ext, proteins_ext) may be empty,
in which case no machine learning or FDR estimation will be performed.
Optional seeds from e.g. untargeted FeatureFinders can be added with
seeds.
Results will be written to features .
Caution: peptide IDs will be shrunk to best hit, FFid metavalues added
and potential seed IDs added.
Usage:
.. code-block:: python
from pyopenms import *
from urllib.request import urlretrieve
urlretrieve("https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/topp/FeatureFinderIdentification_1_input.mzML", "FeatureFinderIdentification_1_input.mzML")
urlretrieve("https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/topp/FeatureFinderIdentification_1_input.idXML", "FeatureFinderIdentification_1_input.idXML")
ffid_algo = FeatureFinderIdentificationAlgorithm()
# load ms data from mzML
mzml = MzMLFile()
mzml_options = mzml.getOptions()
mzml_options.addMSLevel(1) # only MS1
mzml.setOptions(mzml_options)
exp = MSExperiment()
mzml.load("FeatureFinderIdentification_1_input.mzML", exp)
ffid_algo.setMSData(exp)
# annotate mzML file
features = FeatureMap()
features.setPrimaryMSRunPath([b"FeatureFinderIdentification_1_input.idXML"], ffid_algo.getMSData())
peptides = []
proteins = []
peptides_ext = []
proteins_ext = []
IdXMLFile().load("FeatureFinderIdentification_1_input.idXML", proteins, peptides)
#"internal" IDs:
ffid_algo.run(peptides, proteins, peptides_ext, proteins_ext, features)
# Terminal output:
# Summary statistics (counting distinct peptides including PTMs):
# 22 peptides identified (22 internal, 0 additional external)
# 16 peptides with features (16 internal, 0 external)
# 6 peptides without features (6 internal, 0 external)
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::FeatureFinderIdentificationAlgorithm& self, OpenMS::PeptideIdentificationList peptides, const std::vector<OpenMS::ProteinIdentification>& proteins, OpenMS::PeptideIdentificationList peptides_ext, std::vector<OpenMS::ProteinIdentification> proteins_ext, OpenMS::FeatureMap& features, const OpenMS::FeatureMap& seeds, const OpenMS::String& spectra_file) { return self.run(peptides, proteins, peptides_ext, proteins_ext, features, seeds, spectra_file); }, "peptides"_a, "proteins"_a, "peptides_ext"_a, "proteins_ext"_a, "features"_a, "seeds"_a, "spectra_file"_a = "", 
            R"doc(
Run feature detection
:param peptides: Vector of identified peptides
:param proteins: Vector of identified proteins
:param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
:param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
:param features: Feature detection results will be added here
:param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders
)doc")
        .def("runOnCandidates", [](OpenMS::FeatureFinderIdentificationAlgorithm& self, OpenMS::FeatureMap& features) { return self.runOnCandidates(features); }, "features"_a, "Run feature detection on identified features (e.g. loaded from an IdXML file)")
        .def("getMSData", [](OpenMS::FeatureFinderIdentificationAlgorithm& self) -> OpenMS::MSExperiment & { return self.getMSData(); }, nb::rv_policy::reference_internal, "Returns ms data as MSExperiment")
        .def("setMSData", [](OpenMS::FeatureFinderIdentificationAlgorithm& self, const OpenMS::MSExperiment& ms_data) { return self.setMSData(ms_data); }, "ms_data"_a, "Sets ms data")
        .def("setMSData", [](OpenMS::FeatureFinderIdentificationAlgorithm& self, OpenMS::MSExperiment& ms_data) { return self.setMSData(ms_data); }, "ms_data"_a, "Sets ms data")
        .def("getChromatograms", [](OpenMS::FeatureFinderIdentificationAlgorithm& self) -> OpenMS::MSExperiment & { return self.getChromatograms(); }, nb::rv_policy::reference_internal, "Returns chromatogram data as MSExperiment")
        .def("getLibrary", [](OpenMS::FeatureFinderIdentificationAlgorithm& self) -> OpenMS::TargetedExperiment & { return self.getLibrary(); }, nb::rv_policy::reference_internal, "Returns constructed assay library")
        ;

    // -----------------------------------------------------------------------
    // FeatureGroupingAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureGroupingAlgorithm, OpenMS::DefaultParamHandler>(m, "FeatureGroupingAlgorithm", 
        R"doc(
Base class for all feature grouping algorithms
DefaultParamHandler
)doc")
        .def("transferSubelements", [](const OpenMS::FeatureGroupingAlgorithm& self, const std::vector<OpenMS::ConsensusMap>& maps, OpenMS::ConsensusMap& out) { return self.transferSubelements(maps, out); }, "maps"_a, "out"_a, "Transfers subelements (grouped features) from input consensus maps to the result consensus map")
        ;

    // -----------------------------------------------------------------------
    // FeatureGroupingAlgorithmLabeled
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureGroupingAlgorithmLabeled, OpenMS::FeatureGroupingAlgorithm>(m, "FeatureGroupingAlgorithmLabeled", 
        R"doc(
A map feature grouping algorithm for labeling techniques with two
labels
FeatureGroupingAlgorithm
)doc")
        .def(nb::init<>())
        .def("group", [](OpenMS::FeatureGroupingAlgorithmLabeled& self, const std::vector<OpenMS::FeatureMap>& maps, OpenMS::ConsensusMap& out) { return self.group(maps, out); }, "maps"_a, "out"_a)
        .def("transferSubelements", [](const OpenMS::FeatureGroupingAlgorithmLabeled& self, const std::vector<OpenMS::ConsensusMap>& maps, OpenMS::ConsensusMap& out) { return self.transferSubelements(maps, out); }, "maps"_a, "out"_a, "Transfers subelements (grouped features) from input consensus maps to the result consensus map")
        ;

    // -----------------------------------------------------------------------
    // FeatureGroupingAlgorithmQT
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureGroupingAlgorithmQT, OpenMS::FeatureGroupingAlgorithm>(m, "FeatureGroupingAlgorithmQT", 
        R"doc(
A feature grouping algorithm for unlabeled data
FeatureGroupingAlgorithm
)doc")
        .def(nb::init<>())
        .def("group", [](OpenMS::FeatureGroupingAlgorithmQT& self, const std::vector<OpenMS::FeatureMap>& maps, OpenMS::ConsensusMap& out) { return self.group(maps, out); }, "maps"_a, "out"_a)
        .def("group", [](OpenMS::FeatureGroupingAlgorithmQT& self, const std::vector<OpenMS::ConsensusMap>& maps, OpenMS::ConsensusMap& out) { return self.group(maps, out); }, "maps"_a, "out"_a)
        .def("transferSubelements", [](const OpenMS::FeatureGroupingAlgorithmQT& self, const std::vector<OpenMS::ConsensusMap>& maps, OpenMS::ConsensusMap& out) { return self.transferSubelements(maps, out); }, "maps"_a, "out"_a, "Transfers subelements (grouped features) from input consensus maps to the result consensus map")
        ;

    // -----------------------------------------------------------------------
    // FeatureGroupingAlgorithmUnlabeled
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureGroupingAlgorithmUnlabeled, OpenMS::FeatureGroupingAlgorithm>(m, "FeatureGroupingAlgorithmUnlabeled", 
        R"doc(
A map feature grouping algorithm for unlabeled data
FeatureGroupingAlgorithm
)doc")
        .def(nb::init<>())
        .def("getResultMap", [](OpenMS::FeatureGroupingAlgorithmUnlabeled& self) -> OpenMS::ConsensusMap & { return self.getResultMap(); }, nb::rv_policy::reference_internal)
        .def("group", [](OpenMS::FeatureGroupingAlgorithmUnlabeled& self, const std::vector<OpenMS::FeatureMap>& maps, OpenMS::ConsensusMap& out) { return self.group(maps, out); }, "maps"_a, "out"_a)
        .def("addToGroup", [](OpenMS::FeatureGroupingAlgorithmUnlabeled& self, int map_id, const OpenMS::FeatureMap& feature_map) { return self.addToGroup(map_id, feature_map); }, "map_id"_a, "feature_map"_a)
        .def("transferSubelements", [](const OpenMS::FeatureGroupingAlgorithmUnlabeled& self, const std::vector<OpenMS::ConsensusMap>& maps, OpenMS::ConsensusMap& out) { return self.transferSubelements(maps, out); }, "maps"_a, "out"_a, "Transfers subelements (grouped features) from input consensus maps to the result consensus map")
        .def("setReference", [](OpenMS::FeatureGroupingAlgorithmUnlabeled& self, int map_id, const OpenMS::FeatureMap& map) { self.setReference(map_id, map); }, "map_id"_a, "map"_a)
        ;

    // -----------------------------------------------------------------------
    // File
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::File>(m, "File", "Basic file handling operations")
        .def_static("getExecutablePath", []() { return OpenMS::File::getExecutablePath(); })
        .def_static("empty", [](const OpenMS::String& file) { return OpenMS::File::empty(file); }, "file"_a)
        .def_static("rename", [](const OpenMS::String& from, const OpenMS::String& to, bool overwrite_existing, bool verbose) { return OpenMS::File::rename(from, to, overwrite_existing, verbose); }, "from"_a, "to"_a, "overwrite_existing"_a, "verbose"_a)
        .def_static("remove", [](const OpenMS::String& file) { return OpenMS::File::remove(file); }, "file"_a)
        .def_static("removeDirRecursively", [](const OpenMS::String& dir_name) { return OpenMS::File::removeDirRecursively(dir_name); }, "dir_name"_a)
        .def_static("find", [](const OpenMS::String& filename, std::vector<OpenMS::String> directories) { return OpenMS::File::find(filename, directories); }, "filename"_a, "directories"_a)
        .def_static("fileList", [](const OpenMS::String& dir, const OpenMS::String& file_pattern, bool full_path) { std::vector<OpenMS::String> output; OpenMS::File::fileList(dir, file_pattern, output, full_path); return output; }, "dir"_a, "file_pattern"_a, "full_path"_a = false, "Returns list of files matching @p file_pattern in @p dir (returns filenames without paths unless @p full_path is true)")
        .def_static("findDoc", [](const OpenMS::String& filename) { return OpenMS::File::findDoc(filename); }, "filename"_a)
        .def_static("getOpenMSDataPath", []() { return OpenMS::File::getOpenMSDataPath(); })
        .def_static("getOpenMSHomePath", []() { return OpenMS::File::getOpenMSHomePath(); })
        .def_static("getSystemParameters", []() { return OpenMS::File::getSystemParameters(); })
        .def_static("findDatabase", [](const OpenMS::String& db_name) { return OpenMS::File::findDatabase(db_name); }, "db_name"_a)
        .def_static("findExecutable", [](OpenMS::String& exe_filename) { return OpenMS::File::findExecutable(exe_filename); }, "exe_filename"_a)
        .def_static("getTemporaryFile", [](const OpenMS::String& alternative_file) { return OpenMS::File::getTemporaryFile(alternative_file); }, "alternative_file"_a)

        .def_static("exists", [](const OpenMS::String& file) {
            return OpenMS::File::exists(file);
        }, "file"_a, "Check if a file exists")

        .def_static("readable", [](const OpenMS::String& file) {
            return OpenMS::File::readable(file);
        }, "file"_a, "Check if a file is readable")

        .def_static("writable", [](const OpenMS::String& file) {
            return OpenMS::File::writable(file);
        }, "file"_a, "Check if a file is writable")

        .def_static("isDirectory", [](const OpenMS::String& path) {
            return OpenMS::File::isDirectory(path);
        }, "path"_a, "Check if a path is a directory")

        .def_static("basename", [](const OpenMS::String& file) {
            return OpenMS::File::basename(file);
        }, "file"_a, "Get the basename of a file path")

        .def_static("path", [](const OpenMS::String& file) {
            return OpenMS::File::path(file);
        }, "file"_a, "Get the directory path of a file")

        .def_static("absolutePath", [](const OpenMS::String& file) {
            return OpenMS::File::absolutePath(file);
        }, "file"_a, "Get the absolute path")

        .def_static("getTempDirectory", []() {
            return OpenMS::File::getTempDirectory();
        }, "Get the temp directory")

        .def_static("getUserDirectory", []() {
            return OpenMS::File::getUserDirectory();
        }, "Get the user home directory")

        .def_static("getUniqueName", [](bool include_hostname) {
            return OpenMS::File::getUniqueName(include_hostname);
        }, "include_hostname"_a = true, "Get a unique name")
        .def_static("getUniqueName", []() {
            return OpenMS::File::getUniqueName();
        }, "Get a unique name")
        ;

    // -----------------------------------------------------------------------
    // Fitter1D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Fitter1D, OpenMS::DefaultParamHandler>(m, "Fitter1D", 
        R"doc(
Abstract base class for all 1D-dimensional model fitter
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Fitter1D &>())
        .def("__copy__", [](const OpenMS::Fitter1D& self) { return OpenMS::Fitter1D(self); })
        .def("__deepcopy__", [](const OpenMS::Fitter1D& self, nb::dict) { return OpenMS::Fitter1D(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // IDDecoyProbability
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDDecoyProbability, OpenMS::DefaultParamHandler>(m, "IDDecoyProbability", 
        R"doc(
IDDecoyProbability calculates probabilities using decoy approach
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IDDecoyProbability &>())
        .def("__copy__", [](const OpenMS::IDDecoyProbability& self) { return OpenMS::IDDecoyProbability(self); })
        .def("__deepcopy__", [](const OpenMS::IDDecoyProbability& self, nb::dict) { return OpenMS::IDDecoyProbability(self); }, "memo"_a)
        .def("apply", [](OpenMS::IDDecoyProbability& self, OpenMS::PeptideIdentificationList& prob_ids, const OpenMS::PeptideIdentificationList& fwd_ids, const OpenMS::PeptideIdentificationList& rev_ids) { return self.apply(prob_ids, fwd_ids, rev_ids); }, "prob_ids"_a, "fwd_ids"_a, "rev_ids"_a)
        .def("apply", [](OpenMS::IDDecoyProbability& self, OpenMS::PeptideIdentificationList& ids) { return self.apply(ids); }, "ids"_a, 
            R"doc(
Converts the forward and reverse identification into probabilities
:param prob_ids: Output of the algorithm which includes identifications with probability based scores
:param fwd_ids: Input parameter which represents the identifications of the forward search
:param rev_ids: Input parameter which represents the identifications of the reversed search
)doc")
        ;

    // -----------------------------------------------------------------------
    // IDMapper
    // -----------------------------------------------------------------------
    auto idmapper_class = nb::class_<OpenMS::IDMapper, OpenMS::DefaultParamHandler>(m, "IDMapper", 
        R"doc(
Annotates an MSExperiment, FeatureMap or ConsensusMap with peptide
identifications
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IDMapper &>())
        .def("__copy__", [](const OpenMS::IDMapper& self) { return OpenMS::IDMapper(self); })
        .def("__deepcopy__", [](const OpenMS::IDMapper& self, nb::dict) { return OpenMS::IDMapper(self); }, "memo"_a)
        .def("annotate", [](OpenMS::IDMapper& self, OpenMS::AnnotatedMSRun& map, const OpenMS::PeptideIdentificationList& peptide_ids, const std::vector<OpenMS::ProteinIdentification>& protein_ids, bool clear_ids, bool map_ms1) { return self.annotate(map, peptide_ids, protein_ids, clear_ids, map_ms1); }, "map"_a, "peptide_ids"_a, "protein_ids"_a, "clear_ids"_a = false, "map_ms1"_a = false)
        .def("annotate", [](OpenMS::IDMapper& self, OpenMS::AnnotatedMSRun& map, const OpenMS::FeatureMap& fmap, bool clear_ids, bool map_ms1) { return self.annotate(map, fmap, clear_ids, map_ms1); }, "map"_a, "fmap"_a, "clear_ids"_a = false, "map_ms1"_a = false, 
            R"doc(
Mapping method for peak maps\n
The identifications stored in a PeptideIdentification instance can be added to the
corresponding spectrum
Note that a PeptideIdentication is added to ALL spectra which are within the allowed RT and MZ boundaries
:param map: AnnotatedMSRun to receive the identifications
:param peptide_ids: PeptideIdentification for the MSExperiment
:param protein_ids: ProteinIdentification for the MSExperiment
:param clear_ids: Reset peptide and protein identifications of each scan before annotating
:param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
:raises:
Exception: MissingInformation is thrown if entries of 'peptide_ids' do not contain 'MZ' and 'RT' information
)doc")
        .def("annotate", [](OpenMS::IDMapper& self, OpenMS::FeatureMap& map, const OpenMS::PeptideIdentificationList& ids, const std::vector<OpenMS::ProteinIdentification>& protein_ids, bool use_centroid_rt, bool use_centroid_mz, const OpenMS::MSExperiment& spectra) { return self.annotate(map, ids, protein_ids, use_centroid_rt, use_centroid_mz, spectra); }, "map"_a, "ids"_a, "protein_ids"_a, "use_centroid_rt"_a = false, "use_centroid_mz"_a = false, "spectra"_a, 
            R"doc(
Mapping method for peak maps\n
Add peptide identifications stored in a feature map to their
corresponding spectrum
This function converts the feature map to a vector of peptide identifications (all peptide IDs from each feature are taken)
and calls the respective annotate() function
RT and m/z are taken from the peptides, or (if missing) from the feature itself
:param map: AnnotatedMSRun to receive the identifications
:param fmap: FeatureMap with PeptideIdentifications for the MSExperiment
:param clear_ids: Reset peptide and protein identifications of each scan before annotating
:param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
)doc")
        .def("annotate", [](OpenMS::IDMapper& self, OpenMS::ConsensusMap& map, const OpenMS::PeptideIdentificationList& ids, const std::vector<OpenMS::ProteinIdentification>& protein_ids, bool measure_from_subelements, bool annotate_ids_with_subelements, const OpenMS::MSExperiment& spectra) { return self.annotate(map, ids, protein_ids, measure_from_subelements, annotate_ids_with_subelements, spectra); }, "map"_a, "ids"_a, "protein_ids"_a, "measure_from_subelements"_a = false, "annotate_ids_with_subelements"_a = false, "spectra"_a, 
            R"doc(
Mapping method for peak maps\n
Add peptide identifications stored in a feature map to their
corresponding spectrum
This function converts the feature map to a vector of peptide identifications (all peptide IDs from each feature are taken)
and calls the respective annotate() function
RT and m/z are taken from the peptides, or (if missing) from the feature itself
:param map: AnnotatedMSRun to receive the identifications
:param fmap: FeatureMap with PeptideIdentifications for the MSExperiment
:param clear_ids: Reset peptide and protein identifications of each scan before annotating
:param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
)doc")
        .def_static("mapPrecursorsToIdentifications", [](const OpenMS::MSExperiment& spectra, const OpenMS::PeptideIdentificationList& ids, double mz_tol, double rt_tol) { return OpenMS::IDMapper::mapPrecursorsToIdentifications(spectra, ids, mz_tol, rt_tol); }, "spectra"_a, "ids"_a, "mz_tol"_a, "rt_tol"_a, 
            R"doc(
Mapping method for peak maps\n
If all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls by default. If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes) if 'use_centroid_rt' or 'use_centroid_mz' are true\n
In any case, tolerance in RT and m/z dimension is applied according to the global parameters 'rt_tolerance' and 'mz_tolerance'. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value\n
If several features (incl. tolerance) overlap the position of a peptide identification, the identification is annotated to all of them
:param map: MSExperiment to receive the identifications
:param ids: PeptideIdentification for the MSExperiment
:param protein_ids: ProteinIdentification for the MSExperiment
:param measure_from_subelements: Boolean operator set to true if distance estimate from FeatureHandles instead of Centroid
:param annotate_ids_with_subelements: Boolean operator set to true if store map index of FeatureHandle in peptide identification
:param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
:raises:
Exception: MissingInformation is thrown if entries of 'ids' do not contain 'MZ' and 'RT' information
)doc")
        ;
    // Measure enum nested under IDMapper
    nb::enum_<OpenMS::IDMapper::Measure>(idmapper_class, "Measure", nb::is_arithmetic(), "Mass tolerance measure (PPM or Dalton)")
        .value("MEASURE_PPM", OpenMS::IDMapper::Measure::MEASURE_PPM)
        .value("MEASURE_DA", OpenMS::IDMapper::Measure::MEASURE_DA)

        .export_values();

    // -----------------------------------------------------------------------
    // RipFileIdentifier (IDRipper::RipFileIdentifier)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDRipper::RipFileIdentifier>(m, "RipFileIdentifier",
        "Identifies an IDRipper output file")
        .def("getIdentRunIdx", &OpenMS::IDRipper::RipFileIdentifier::getIdentRunIdx, "Get identification run index")
        .def("getFileOriginIdx", &OpenMS::IDRipper::RipFileIdentifier::getFileOriginIdx, "Get file origin index")
        .def("getOriginFullname", &OpenMS::IDRipper::RipFileIdentifier::getOriginFullname, nb::rv_policy::reference_internal, "Get origin full name")
        .def("getOutputBasename", &OpenMS::IDRipper::RipFileIdentifier::getOutputBasename, nb::rv_policy::reference_internal, "Get output base name")
        ;

    // -----------------------------------------------------------------------
    // RipFileContent (IDRipper::RipFileContent)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDRipper::RipFileContent>(m, "RipFileContent",
        "Represents the content of an IDRipper output file")
        .def("getProteinIdentifications", &OpenMS::IDRipper::RipFileContent::getProteinIdentifications, nb::rv_policy::reference_internal, "Get protein identifications")
        .def("getPeptideIdentifications", &OpenMS::IDRipper::RipFileContent::getPeptideIdentifications, nb::rv_policy::reference_internal, "Get peptide identifications")
        ;

    // -----------------------------------------------------------------------
    // IDRipper
    // -----------------------------------------------------------------------
    auto idripper_class = nb::class_<OpenMS::IDRipper, OpenMS::DefaultParamHandler>(m, "IDRipper", 
        R"doc(
Ripping protein/peptide identification according their file origin
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("rip", [](OpenMS::IDRipper& self, std::vector<OpenMS::ProteinIdentification> proteins, OpenMS::PeptideIdentificationList& peptides, bool numeric_filenames, bool split_ident_runs) { std::map<OpenMS::IDRipper::RipFileIdentifier, OpenMS::IDRipper::RipFileContent, OpenMS::IDRipper::RipFileIdentifierIdxComparator> ripped; self.rip(ripped, proteins, peptides, numeric_filenames, split_ident_runs); return ripped; }, "proteins"_a, "peptides"_a, "numeric_filenames"_a, "split_ident_runs"_a)
        .def("rip", [](OpenMS::IDRipper& self, std::vector<OpenMS::ProteinIdentification> proteins, OpenMS::PeptideIdentificationList& peptides, bool numeric_filenames, bool split_ident_runs, bool /*use_vector_output*/) { std::vector<OpenMS::IDRipper::RipFileIdentifier> rfis; std::vector<OpenMS::IDRipper::RipFileContent> rfcs; self.rip(rfis, rfcs, proteins, peptides, numeric_filenames, split_ident_runs); return nb::make_tuple(rfis, rfcs); }, "proteins"_a, "peptides"_a, "numeric_filenames"_a, "split_ident_runs"_a, "use_vector_output"_a)
        ;
    // OriginAnnotationFormat enum nested under IDRipper
    nb::enum_<OpenMS::IDRipper::OriginAnnotationFormat>(idripper_class, "OriginAnnotationFormat", nb::is_arithmetic(), "Format for annotating the origin of identifications")
        .value("FILE_ORIGIN", OpenMS::IDRipper::OriginAnnotationFormat::FILE_ORIGIN)
        .value("MAP_INDEX", OpenMS::IDRipper::OriginAnnotationFormat::MAP_INDEX)
        .value("ID_MERGE_INDEX", OpenMS::IDRipper::OriginAnnotationFormat::ID_MERGE_INDEX)
        .value("UNKNOWN_OAF", OpenMS::IDRipper::OriginAnnotationFormat::UNKNOWN_OAF)
        .value("SIZE_OF_ORIGIN_ANNOTATION_FORMAT", OpenMS::IDRipper::OriginAnnotationFormat::SIZE_OF_ORIGIN_ANNOTATION_FORMAT)

        .export_values();

    // -----------------------------------------------------------------------
    // IDScoreSwitcherAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDScoreSwitcherAlgorithm, OpenMS::DefaultParamHandler>(m, "IDScoreSwitcherAlgorithm", 
        R"doc(
DefaultParamHandler

Algorithm to switch identification scores within identification or consensus feature maps
This class provides functionality to switch the main scoring type used in peptide or protein
identification data. It supports switching between different score types, such as raw scores,
E-values, posterior probabilities, posterior error probabilities, FDR, and q-values.
)doc")
        .def(nb::init<>())
        .def("isScoreType", [](const OpenMS::IDScoreSwitcherAlgorithm& self, const OpenMS::String& score_name, const OpenMS::Scores::IDType& type) { return self.isScoreType(score_name, type); }, "score_name"_a, "type"_a)
        .def_static("toScoreTypeEnum", [](const OpenMS::String& score_type) { return OpenMS::IDScoreSwitcherAlgorithm::toScoreTypeEnum(score_type); }, "score_type"_a, 
            R"doc(
Searches for a score type in a PeptideIdentification
Returns a ScoreSearchResult indicating whether the main score is of the
requested type, and if not, searches for scores of that type in the
meta values of the first hit.
:param id: The PeptideIdentification to analyze
:param score_type: The IDType to search for (e.g., IDType.PEP)
:returns: ScoreSearchResult with is_main_score_type and score_name fields
)doc")
        .def("isScoreTypeHigherBetter", [](OpenMS::IDScoreSwitcherAlgorithm& self, OpenMS::Scores::IDType score_type) { return self.isScoreTypeHigherBetter(score_type); }, "score_type"_a, 
            R"doc(
Converts a string representation of a score type to an IDType enum
:param score_type: The string representation of the score type
:returns: The corresponding IDType enum value
:raises: Exception::MissingInformation if the score_type string is not recognized
)doc")
        .def("getScoreNames", [](OpenMS::IDScoreSwitcherAlgorithm& self) { return self.getScoreNames(); }, 
            R"doc(
Determines whether a higher score type is better given an IDType enum
:param score_type: The score type to check
:returns: True if a higher score type is better
)doc")
        .def("switchToGeneralScoreType", [](OpenMS::IDScoreSwitcherAlgorithm& self, OpenMS::PeptideIdentificationList& pep_ids, OpenMS::Scores::IDType type, size_t counter) { self.switchToGeneralScoreType(pep_ids, type, counter); return counter; }, "pep_ids"_a, "type"_a, "counter"_a,
            R"doc(
Gets a vector of all score names that are used in OpenMS
:returns: A vector of all score names (e.g., "q-value", "ln(hyperscore)")
)doc")
        .def("switchToGeneralScoreType", [](OpenMS::IDScoreSwitcherAlgorithm& self, OpenMS::ConsensusMap& cmap, OpenMS::Scores::IDType type, size_t counter, bool unassigned_peptides_too) { self.switchToGeneralScoreType(cmap, type, counter, unassigned_peptides_too); return counter; }, "cmap"_a, "type"_a, "counter"_a, "unassigned_peptides_too"_a = true,
            R"doc(
Switches the score type of a PeptideIdentificationList to a general score type
:param pep_ids: The PeptideIdentificationList whose scores need to be switched
:param type: The desired general score type to switch to
:param counter: A reference to a counter that will be incremented for each peptide identification processed
)doc")
        .def("switchScores", [](OpenMS::IDScoreSwitcherAlgorithm& self, OpenMS::ConsensusMap& cmap, size_t counter, bool unassigned_peptides_too) { self.switchScores(cmap, counter, unassigned_peptides_too); return counter; }, "cmap"_a, "counter"_a, "unassigned_peptides_too"_a = true,
            R"doc(
Switches the scores of peptide identifications
:param pep_ids: The peptide identifications whose scores need to be switched
:param counter: A reference to a counter that will be incremented for each peptide identification processed
)doc")
        .def("switchScores", [](OpenMS::IDScoreSwitcherAlgorithm& self, OpenMS::PeptideIdentificationList& pep_ids, size_t counter) { self.switchScores(pep_ids, counter); return counter; }, "pep_ids"_a, "counter"_a,
            R"doc(
Switches the score type of a ConsensusMap to a general score type
:param cmap: The ConsensusMap containing peptide identifications whose scores need to be switched
:param type: The desired general score type to switch to
:param counter: A reference to a counter that will be incremented for each peptide identification processed
:param unassigned_peptides_too: Whether to include unassigned peptides in the score switching process
)doc")

        .def("findScoreType", [](OpenMS::IDScoreSwitcherAlgorithm& self, OpenMS::PeptideIdentification& id, OpenMS::Scores::IDType score_type) {
            return self.findScoreType(id, score_type);
        }, "id"_a, "score_type"_a, "Search for a score type in a PeptideIdentification")

        .def_static("switchToScoreType", nb::overload_cast<OpenMS::PeptideIdentificationList&, OpenMS::String>(&OpenMS::IDScoreSwitcherAlgorithm::switchToScoreType),
            "pep_ids"_a, "requested_score_type_as_string"_a, "Switch the score type of peptide identifications")
        .def_static("switchToScoreType", [](OpenMS::ConsensusMap& cmap, OpenMS::String requested_score_type, bool include_unassigned) {
            return OpenMS::IDScoreSwitcherAlgorithm::switchToScoreType(cmap, requested_score_type, include_unassigned);
        }, "cmap"_a, "requested_score_type_as_string"_a, "include_unassigned"_a = true, "Switch the score type of a ConsensusMap")

        .def_static("switchBackScoreType", nb::overload_cast<OpenMS::PeptideIdentificationList&, OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult>(&OpenMS::IDScoreSwitcherAlgorithm::switchBackScoreType),
            "pep_ids"_a, "isr"_a, "Revert scores to original type")
        .def_static("switchBackScoreType", [](OpenMS::ConsensusMap& cmap, OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult isr, bool include_unassigned) {
            OpenMS::IDScoreSwitcherAlgorithm::switchBackScoreType(cmap, isr, include_unassigned);
        }, "cmap"_a, "isr"_a, "include_unassigned"_a = true, "Revert scores of a ConsensusMap to original type")
        ;


    // -----------------------------------------------------------------------
    // ScoreSearchResult
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult>(m, "ScoreSearchResult", "OpenMS class ScoreSearchResult")
        .def(nb::init<>())
        .def_rw("is_main_score_type", &OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult::is_main_score_type,
                "True if the main score is already of the requested score type")
        .def_rw("score_name", &OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult::score_name,
                "Name of score to use (main score name or meta value name)")
        .def("__repr__", [](const OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult& self) {
            std::ostringstream oss;
            oss << "ScoreSearchResult(is_main_score_type=" << (self.is_main_score_type ? "True" : "False")
                << ", score_name='" << self.score_name << "')";
            return oss.str();
        })
        ;


    // -----------------------------------------------------------------------
    // IDSwitchResult
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult>(m, "IDSwitchResult", "OpenMS class IDSwitchResult")
        .def(nb::init<>())
        .def_rw("original_score_name", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::original_score_name,
                "The name of the original score used before the switch")
        .def_rw("original_score_higher_better", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::original_score_higher_better,
                "Whether a higher original score is better")
        .def_rw("original_score_type", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::original_score_type,
                "The type of the original score")
        .def_rw("requested_score_higher_better", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::requested_score_higher_better,
                "Whether a higher requested score is better")
        .def_rw("requested_score_type", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::requested_score_type,
                "The type of the requested score")
        .def_rw("requested_score_name", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::requested_score_name,
                "The search engine score name (e.g. 'X!Tandem_score') or score category (e.g. 'PEP')")
        .def_rw("score_switched", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::score_switched,
                "Flag indicating whether the main score was switched")
        .def("__repr__", [](const OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult& self) {
            std::ostringstream oss;
            oss << "IDSwitchResult(original_score_name='" << self.original_score_name
                << "', original_score_higher_better=" << (self.original_score_higher_better ? "True" : "False")
                << ", score_switched=" << (self.score_switched ? "True" : "False")
                << ")";
            return oss.str();
        })
        ;

    // -----------------------------------------------------------------------
    // IMSDataConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Interfaces::IMSDataConsumer>(m, "IMSDataConsumer", "The interface of a consumer of spectra and chromatograms")
        .def("consumeSpectrum", [](OpenMS::Interfaces::IMSDataConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a, "Consume a spectrum. The spectrum will be consumed by the implementation and possibly modified")
        .def("consumeChromatogram", [](OpenMS::Interfaces::IMSDataConsumer& self, OpenMS::MSChromatogram& c) { return self.consumeChromatogram(c); }, "c"_a, "Consume a chromatogram. The chromatogram will be consumed by the implementation and possibly modified")
        .def("setExpectedSize", [](OpenMS::Interfaces::IMSDataConsumer& self, size_t expectedSpectra, size_t expectedChromatograms) { return self.setExpectedSize(expectedSpectra, expectedChromatograms); }, "expectedSpectra"_a, "expectedChromatograms"_a)
        .def("setExperimentalSettings", [](OpenMS::Interfaces::IMSDataConsumer& self, const OpenMS::ExperimentalSettings& exp) { return self.setExperimentalSettings(exp); }, "exp"_a, 
            R"doc(
Set expected size of spectra and chromatograms to be consumed\n
Some implementations might care about the number of spectra and
chromatograms to be consumed and need to be informed about this
(usually before consuming starts)
:param expectedSpectra: Number of spectra expected
:param expectedChromatograms: Number of chromatograms expected
)doc")
        ;

    // -----------------------------------------------------------------------
    // IMTypes
    // -----------------------------------------------------------------------
    auto imtypes_class = nb::class_<OpenMS::IMTypes>(m, "IMTypes", "OpenMS class IMTypes")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IMTypes &>())
        .def("__copy__", [](const OpenMS::IMTypes& self) { return OpenMS::IMTypes(self); })
        .def("__deepcopy__", [](const OpenMS::IMTypes& self, nb::dict) { return OpenMS::IMTypes(self); }, "memo"_a)
        .def_static("determineIMFormat", [](const OpenMS::MSExperiment& exp, int ms_level) { return OpenMS::IMTypes::determineIMFormat(exp, ms_level); }, "exp"_a, "ms_level"_a)
        .def_static("determineIMFormat", [](const OpenMS::MSSpectrum& spec) { return OpenMS::IMTypes::determineIMFormat(spec); }, "spec"_a)

        .def_static("toDriftTimeUnit", [](const OpenMS::String& dtu_string) {
            return OpenMS::toDriftTimeUnit(dtu_string);
        }, "dtu_string"_a, "Convert string to DriftTimeUnit")

        .def_static("driftTimeUnitToString", [](OpenMS::DriftTimeUnit value) {
            return OpenMS::driftTimeUnitToString(value);
        }, "value"_a, "Convert DriftTimeUnit to string")

        .def_static("toIMFormat", [](const OpenMS::String& im_format) {
            return OpenMS::toIMFormat(im_format);
        }, "im_format"_a, "Convert string to IMFormat")

        .def_static("imFormatToString", [](OpenMS::IMFormat value) {
            return OpenMS::imFormatToString(value);
        }, "value"_a, "Convert IMFormat to string")

        .def_static("toIMPeakType", [](const OpenMS::String& im_peak_type) {
            return OpenMS::toIMPeakType(im_peak_type);
        }, "im_peak_type"_a, "Convert string to IMPeakType")

        .def_static("imPeakTypeToString", [](OpenMS::IMPeakType value) {
            return OpenMS::imPeakTypeToString(value);
        }, "value"_a, "Convert IMPeakType to string")
        ;

    // -----------------------------------------------------------------------
    // IsobaricQuantifier
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsobaricQuantifier, OpenMS::DefaultParamHandler>(m, "IsobaricQuantifier", 
        R"doc(
Performs isotope correction and normalization of isobaric labeling data.
Given extracted channel intensities, this class corrects for isotope
impurities using a correction matrix and optionally normalizes the
intensities for further downstream processing
DefaultParamHandler
)doc")
        .def(nb::init<const OpenMS::IsobaricQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::IsobaricQuantifier &>())
        .def("__copy__", [](const OpenMS::IsobaricQuantifier& self) { return OpenMS::IsobaricQuantifier(self); })
        .def("__deepcopy__", [](const OpenMS::IsobaricQuantifier& self, nb::dict) { return OpenMS::IsobaricQuantifier(self); }, "memo"_a)
        .def("quantify", [](OpenMS::IsobaricQuantifier& self, const OpenMS::ConsensusMap& consensus_map_in, OpenMS::ConsensusMap& consensus_map_out) { self.quantify(consensus_map_in, consensus_map_out); }, "consensus_map_in"_a, "consensus_map_out"_a, "Quantifies isobaric labeled peptides/proteins")
        ;

    // -----------------------------------------------------------------------
    // IsobaricQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsobaricQuantitationMethod, OpenMS::DefaultParamHandler>(m, "IsobaricQuantitationMethod", 
        R"doc(
Abstract base class describing an isobaric quantitation method in terms of
the reporter ion channels used and an isotope correction matrix. Isobaric
labeling methods like TMT and iTRAQ use reporter ions for multiplexed
quantitation of peptides/proteins across multiple samples
DefaultParamHandler
)doc")
        .def("getChannelInformation", [](const OpenMS::IsobaricQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::IsobaricQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getReferenceChannel", [](const OpenMS::IsobaricQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")

        .def("getIsotopeCorrectionMatrix", [](const OpenMS::IsobaricQuantitationMethod& self) {
            return self.getIsotopeCorrectionMatrix();
        }, "Get the isotope correction matrix")
        ;

    // -----------------------------------------------------------------------
    // IsotopeLabelingMDVs
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsotopeLabelingMDVs, OpenMS::DefaultParamHandler>(m, "IsotopeLabelingMDVs", 
        R"doc(
IsotopeLabelingMDVs is a class to support and analyze isotopic
labeling experiments (i.e. MDVs : Mass Distribution Vectors, also
known as Mass Isotopomer Distribution (MID))
)doc")
        .def(nb::init<>())
        .def("isotopicCorrection", [](OpenMS::IsotopeLabelingMDVs& self, const OpenMS::Feature& normalized_feature, OpenMS::Feature& corrected_feature, const OpenMS::Matrix<double>& correction_matrix, const OpenMS::IsotopeLabelingMDVs::DerivatizationAgent& correction_matrix_agent) { return self.isotopicCorrection(normalized_feature, corrected_feature, correction_matrix, correction_matrix_agent); }, "normalized_feature"_a, "corrected_feature"_a, "correction_matrix"_a, "correction_matrix_agent"_a)
        .def("isotopicCorrections", [](OpenMS::IsotopeLabelingMDVs& self, const OpenMS::FeatureMap& measured_fm, OpenMS::FeatureMap& corrected_fm, const OpenMS::Matrix<double>& correction_matrix, const OpenMS::IsotopeLabelingMDVs::DerivatizationAgent& correction_matrix_agent) { return self.isotopicCorrections(measured_fm, corrected_fm, correction_matrix, correction_matrix_agent); }, "measured_fm"_a, "corrected_fm"_a, "correction_matrix"_a, "correction_matrix_agent"_a, 
            R"doc(
This function performs an isotopic correction to account for unlabeled abundances coming from
the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
:param normalized_feature: Feature with normalized values for each component and unlabeled chemical formula for each component group
:param correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
:param correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
:return: corrected_feature: Feature with corrected values for each component
)doc")
        .def("calculateIsotopicPurity", [](OpenMS::IsotopeLabelingMDVs& self, const std::vector<double>& experiment_data, const OpenMS::String& isotopic_purity_name) { OpenMS::Feature normalized_feature; self.calculateIsotopicPurity(normalized_feature, experiment_data, isotopic_purity_name); return normalized_feature; }, "experiment_data"_a, "isotopic_purity_name"_a, 
            R"doc(
This function performs an isotopic correction to account for unlabeled abundances coming from
the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
:param normalized_featuremap: FeatureMap with normalized values for each component and unlabeled chemical formula for each component group
:param correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
:param correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
:return corrected_featuremap: FeatureMap with corrected values for each component
)doc")
        .def("calculateMDVAccuracy", [](OpenMS::IsotopeLabelingMDVs& self, const OpenMS::String& feature_name, const OpenMS::String& fragment_isotopomer_theoretical_formula) { OpenMS::Feature normalized_feature; self.calculateMDVAccuracy(normalized_feature, feature_name, fragment_isotopomer_theoretical_formula); return normalized_feature; }, "feature_name"_a, "fragment_isotopomer_theoretical_formula"_a, 
            R"doc(
This function calculates the isotopic purity of the MDV using the following formula:
isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
where n in M+n is represented as the index of the result
The formula is extracted from "High-resolution 13C metabolic flux analysis",
Long et al, doi:10.1038/s41596-019-0204-0
:param normalized_feature: Feature with normalized values for each component and the number of heavy labeled e.g., carbons. Out is a Feature with the calculated isotopic purity for the component group
:param experiment_data: Vector of experiment data in percent
:param isotopic_purity_name: Name of the isotopic purity tracer to be saved as a meta value
)doc")
        .def("calculateMDV", [](OpenMS::IsotopeLabelingMDVs& self, const OpenMS::Feature& measured_feature, const OpenMS::IsotopeLabelingMDVs::MassIntensityType& mass_intensity_type, const OpenMS::String& feature_name) { OpenMS::Feature normalized_feature; self.calculateMDV(measured_feature, normalized_feature, mass_intensity_type, feature_name); return normalized_feature; }, "measured_feature"_a, "mass_intensity_type"_a, "feature_name"_a, 
            R"doc(
This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
using average deviation to the mean
param normalized_featuremap: FeatureMap with normalized values for each component and the chemical formula of the component group. Out is a FeatureMap with the component group accuracy and accuracy for the error for each component
param fragment_isotopomer_measured: Measured scan values
param fragment_isotopomer_theoretical_formula: A map of ProteinName/peptideRef to Empirical formula from which the theoretical values will be generated
)doc")
        .def("calculateMDVs", [](OpenMS::IsotopeLabelingMDVs& self, const OpenMS::FeatureMap& measured_fm, const OpenMS::IsotopeLabelingMDVs::MassIntensityType& mass_intensity_type, const OpenMS::String& feature_name) { OpenMS::FeatureMap normalized_fm; self.calculateMDVs(measured_fm, normalized_fm, mass_intensity_type, feature_name); return normalized_fm; }, "measured_fm"_a, "mass_intensity_type"_a, "feature_name"_a)
        .def("calculateMDVAccuracies", [](OpenMS::IsotopeLabelingMDVs& self, OpenMS::FeatureMap& normalized_fm, const std::string& feature_name, const std::map<std::string, std::string>& fragment_isotopomer_theoretical_formulas) { self.calculateMDVAccuracies(normalized_fm, feature_name, fragment_isotopomer_theoretical_formulas); }, "normalized_fm"_a, "feature_name"_a, "fragment_isotopomer_theoretical_formulas"_a, "Calculate MDV accuracies")
        ;

    // -----------------------------------------------------------------------
    // ItraqEightPlexQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ItraqEightPlexQuantitationMethod, OpenMS::IsobaricQuantitationMethod>(m, "ItraqEightPlexQuantitationMethod", 
        R"doc(
iTRAQ 8 plex quantitation to be used with the IsobaricQuantitation
IsobaricQuantitationMethod
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ItraqEightPlexQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::ItraqEightPlexQuantitationMethod& self) { return OpenMS::ItraqEightPlexQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::ItraqEightPlexQuantitationMethod& self, nb::dict) { return OpenMS::ItraqEightPlexQuantitationMethod(self); }, "memo"_a)
        .def("getChannelInformation", [](const OpenMS::ItraqEightPlexQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::ItraqEightPlexQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::ItraqEightPlexQuantitationMethod& self) { return self.getIsotopeCorrectionMatrix(); }, "Returns the isotope correction matrix for correcting reporter ion intensities")
        .def("getReferenceChannel", [](const OpenMS::ItraqEightPlexQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")
        ;

    // -----------------------------------------------------------------------
    // ItraqFourPlexQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ItraqFourPlexQuantitationMethod, OpenMS::IsobaricQuantitationMethod>(m, "ItraqFourPlexQuantitationMethod", 
        R"doc(
iTRAQ 4 plex quantitation to be used with the IsobaricQuantitation
IsobaricQuantitationMethod
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ItraqFourPlexQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::ItraqFourPlexQuantitationMethod& self) { return OpenMS::ItraqFourPlexQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::ItraqFourPlexQuantitationMethod& self, nb::dict) { return OpenMS::ItraqFourPlexQuantitationMethod(self); }, "memo"_a)
        .def("getChannelInformation", [](const OpenMS::ItraqFourPlexQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::ItraqFourPlexQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::ItraqFourPlexQuantitationMethod& self) { return self.getIsotopeCorrectionMatrix(); }, "Returns the isotope correction matrix for correcting reporter ion intensities")
        .def("getReferenceChannel", [](const OpenMS::ItraqFourPlexQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")
        ;

    // -----------------------------------------------------------------------
    // JavaInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::JavaInfo>(m, "JavaInfo", "Detect Java and retrieve information")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::JavaInfo &>())
        .def("__copy__", [](const OpenMS::JavaInfo& self) { return OpenMS::JavaInfo(self); })
        .def("__deepcopy__", [](const OpenMS::JavaInfo& self, nb::dict) { return OpenMS::JavaInfo(self); }, "memo"_a)
        .def_static("canRun", [](const OpenMS::String& java_executable, bool verbose_on_error) { return OpenMS::JavaInfo::canRun(java_executable, verbose_on_error); }, "java_executable"_a, "verbose_on_error"_a)
        ;

    // -----------------------------------------------------------------------
    // KDTreeFeatureMaps
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::KDTreeFeatureMaps, OpenMS::DefaultParamHandler>(m, "KDTreeFeatureMaps", "OpenMS class KDTreeFeatureMaps")
        .def(nb::init<>())
        .def("rt", [](const OpenMS::KDTreeFeatureMaps& self, size_t i) { return self.rt(i); }, "i"_a)
        .def("mz", [](const OpenMS::KDTreeFeatureMaps& self, size_t i) { return self.mz(i); }, "i"_a)
        .def("intensity", [](const OpenMS::KDTreeFeatureMaps& self, size_t i) { return self.intensity(i); }, "i"_a)
        .def("charge", [](const OpenMS::KDTreeFeatureMaps& self, size_t i) { return self.charge(i); }, "i"_a)
        .def("mapIndex", [](const OpenMS::KDTreeFeatureMaps& self, size_t i) { return self.mapIndex(i); }, "i"_a)
        .def("size", [](const OpenMS::KDTreeFeatureMaps& self) { return self.size(); })
        .def("treeSize", [](const OpenMS::KDTreeFeatureMaps& self) { return self.treeSize(); })
        .def("numMaps", [](const OpenMS::KDTreeFeatureMaps& self) { return self.numMaps(); })
        .def("clear", [](OpenMS::KDTreeFeatureMaps& self) { return self.clear(); })
        .def("optimizeTree", [](OpenMS::KDTreeFeatureMaps& self) { return self.optimizeTree(); })
        .def("getNeighborhood", [](const OpenMS::KDTreeFeatureMaps& self, size_t index, double rt_tol, double mz_tol, bool mz_ppm, bool include_features_from_same_map, double max_pairwise_log_fc) { std::vector<size_t> result_indices; self.getNeighborhood(index, result_indices, rt_tol, mz_tol, mz_ppm, include_features_from_same_map, max_pairwise_log_fc); return result_indices; }, "index"_a, "rt_tol"_a, "mz_tol"_a, "mz_ppm"_a, "include_features_from_same_map"_a, "max_pairwise_log_fc"_a, "Fill `result` with indices of all features compatible (wrt. RT, m/z, map index) to the feature with `index`")
        .def("queryRegion", [](const OpenMS::KDTreeFeatureMaps& self, double rt_low, double rt_high, double mz_low, double mz_high, size_t ignored_map_index) { std::vector<size_t> result_indices; self.queryRegion(rt_low, rt_high, mz_low, mz_high, result_indices, ignored_map_index); return result_indices; }, "rt_low"_a, "rt_high"_a, "mz_low"_a, "mz_high"_a, "ignored_map_index"_a)
        .def("__len__", [](OpenMS::KDTreeFeatureMaps& self) { return self.size(); })
        .def("addMaps", [](OpenMS::KDTreeFeatureMaps& self, const std::vector<OpenMS::FeatureMap>& maps) { self.addMaps(maps); }, "maps"_a)
        .def("applyTransformations", [](OpenMS::KDTreeFeatureMaps& self, const std::vector<OpenMS::TransformationModelLowess*>& trafos) { self.applyTransformations(trafos); }, "trafos"_a, "Apply transformations to RT values")
        ;

    // -----------------------------------------------------------------------
    // LevMarqFitter1D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::LevMarqFitter1D, OpenMS::Fitter1D>(m, "LevMarqFitter1D", "Abstract class for 1D-model fitter using Levenberg-Marquardt algorithm for parameter optimization")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::LevMarqFitter1D &>())
        .def("__copy__", [](const OpenMS::LevMarqFitter1D& self) { return OpenMS::LevMarqFitter1D(self); })
        .def("__deepcopy__", [](const OpenMS::LevMarqFitter1D& self, nb::dict) { return OpenMS::LevMarqFitter1D(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // EmgFitter1D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::EmgFitter1D, OpenMS::LevMarqFitter1D>(m, "EmgFitter1D", 
        R"doc(
Exponentially modified gaussian distribution fitter (1-dim.) using
Levenberg-Marquardt algorithm (Eigen implementation) for parameter
optimization
LevMarqFitter1D
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::EmgFitter1D &>())
        .def("__copy__", [](const OpenMS::EmgFitter1D& self) { return OpenMS::EmgFitter1D(self); })
        .def("__deepcopy__", [](const OpenMS::EmgFitter1D& self, nb::dict) { return OpenMS::EmgFitter1D(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // LowessSmoothing
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::LowessSmoothing, OpenMS::DefaultParamHandler>(m, "LowessSmoothing", 
        R"doc(
LOWESS (locally weighted scatterplot smoothing)
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("smoothData", [](OpenMS::LowessSmoothing& self, const std::vector<double>& x_vals, const std::vector<double>& y_vals) { std::vector<double> smoothed; self.smoothData(x_vals, y_vals, smoothed); return smoothed; }, "x_vals"_a, "y_vals"_a, "Smoothing method that receives x and y coordinates (e.g., RT and intensities) and returns smoothed intensities")
        ;

    // -----------------------------------------------------------------------
    // MRMFeatureFilter
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureFilter, OpenMS::DefaultParamHandler>(m, "MRMFeatureFilter", 
        R"doc(
DefaultParamHandler

Flags or filters MRM features that do not pass QC criteria
This class provides comprehensive quality control filtering for MRM/SRM features.
It can filter based on:
- Retention time, intensity, and quality bounds
- Ion ratios between transitions
- %RSD from pooled QC samples
- Background interference from blank samples
Example usage:
.. code-block:: python
from pyopenms import *
# Load features and transitions
features = FeatureMap()
FeatureXMLFile().load("features.featureXML", features)
transitions = TargetedExperiment()
TraMLFile().load("transitions.traML", transitions)
# Set up QC criteria
qc = MRMFeatureQC()
comp_qc = MRMFeatureQC.ComponentQCs()
comp_qc.component_name = "my_transition"
comp_qc.retention_time_l = 5.0
comp_qc.retention_time_u = 7.0
comp_qc.intensity_l = 1000.0
qc.component_qcs.append(comp_qc)
# Apply filtering
filter = MRMFeatureFilter()
filter.FilterFeatureMap(features, qc, transitions)
)doc")
        .def(nb::init<>())
        .def("FilterFeatureMap", [](OpenMS::MRMFeatureFilter& self, OpenMS::FeatureMap& features, const OpenMS::MRMFeatureQC& filter_criteria, const OpenMS::TargetedExperiment& transitions) { return self.FilterFeatureMap(features, filter_criteria, transitions); }, "features"_a, "filter_criteria"_a, "transitions"_a)
        .def("FilterFeatureMapPercRSD", [](OpenMS::MRMFeatureFilter& self, OpenMS::FeatureMap& features, const OpenMS::MRMFeatureQC& filter_criteria, const OpenMS::MRMFeatureQC& filter_values) { return self.FilterFeatureMapPercRSD(features, filter_criteria, filter_values); }, "features"_a, "filter_criteria"_a, "filter_values"_a, 
            R"doc(
Flags or filters features and subordinates in a FeatureMap
:param features: FeatureMap to flag or filter
:param filter_criteria: MRMFeatureQC class defining QC parameters
:param transitions: Transitions from a TargetedExperiment
)doc")
        .def("FilterFeatureMapBackgroundInterference", [](OpenMS::MRMFeatureFilter& self, OpenMS::FeatureMap& features, const OpenMS::MRMFeatureQC& filter_criteria, const OpenMS::MRMFeatureQC& filter_values) { return self.FilterFeatureMapBackgroundInterference(features, filter_criteria, filter_values); }, "features"_a, "filter_criteria"_a, "filter_values"_a, 
            R"doc(
Filter features based on percent relative standard deviation (%RSD)
Uses QC samples to calculate %RSD and filters features exceeding thresholds.
:param features: FeatureMap to filter
:param filter_criteria: MRMFeatureQC with %RSD thresholds
:param filter_values: MRMFeatureQC with calculated %RSD values
)doc")
        .def("calculateIonRatio", [](const OpenMS::MRMFeatureFilter& self, const OpenMS::Feature& component_1, const OpenMS::Feature& component_2, const OpenMS::String& feature_name) { return self.calculateIonRatio(component_1, component_2, feature_name); }, "component_1"_a, "component_2"_a, "feature_name"_a, 
            R"doc(
Estimate background interference from blank samples
Analyzes blank samples to determine background signal levels.
:param samples: Vector of FeatureMaps from blank samples
:param filter_template: MRMFeatureQC to populate with interference values
:param transitions: TargetedExperiment with transition information
)doc")
        .def("calculateRTDifference", [](const OpenMS::MRMFeatureFilter& self, OpenMS::Feature& component_1, OpenMS::Feature& component_2) { return self.calculateRTDifference(component_1, component_2); }, "component_1"_a, "component_2"_a, 
            R"doc(
Calculate the ion ratio between two features
:param component_1: First feature (e.g., quantifier transition)
:param component_2: Second feature (e.g., qualifier transition)
:param feature_name: Name of the feature to use for ratio (e.g., "peak_apex_int")
:returns: The ratio of component_1 / component_2
)doc")
        .def("calculateResolution", [](const OpenMS::MRMFeatureFilter& self, OpenMS::Feature& component_1, OpenMS::Feature& component_2) { return self.calculateResolution(component_1, component_2); }, "component_1"_a, "component_2"_a, 
            R"doc(
Calculate the retention time difference between two features
:param component_1: First feature
:param component_2: Second feature
:returns: The absolute RT difference in seconds
)doc")
        .def("EstimateDefaultMRMFeatureQCValues", [](const OpenMS::MRMFeatureFilter& self, const std::vector<OpenMS::FeatureMap>& samples, OpenMS::MRMFeatureQC& filter_template, const OpenMS::TargetedExperiment& transitions, const bool& init_template_values) { self.EstimateDefaultMRMFeatureQCValues(samples, filter_template, transitions, init_template_values); }, "samples"_a, "filter_template"_a, "transitions"_a, "init_template_values"_a)
        .def("TransferLLOQAndULOQToCalculatedConcentrationBounds", [](OpenMS::MRMFeatureFilter& self, const std::vector<OpenMS::AbsoluteQuantitationMethod>& quantitation_method, OpenMS::MRMFeatureQC& filter_template) { self.TransferLLOQAndULOQToCalculatedConcentrationBounds(quantitation_method, filter_template); }, "quantitation_method"_a, "filter_template"_a)
        .def("EstimatePercRSD", [](const OpenMS::MRMFeatureFilter& self, const std::vector<OpenMS::FeatureMap>& samples, OpenMS::MRMFeatureQC& filter_template, const OpenMS::TargetedExperiment& transitions) { self.EstimatePercRSD(samples, filter_template, transitions); }, "samples"_a, "filter_template"_a, "transitions"_a)
        .def("EstimateBackgroundInterferences", [](const OpenMS::MRMFeatureFilter& self, const std::vector<OpenMS::FeatureMap>& samples, OpenMS::MRMFeatureQC& filter_template, const OpenMS::TargetedExperiment& transitions) { self.EstimateBackgroundInterferences(samples, filter_template, transitions); }, "samples"_a, "filter_template"_a, "transitions"_a)
        ;

    // -----------------------------------------------------------------------
    // MRMMapping
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMMapping, OpenMS::DefaultParamHandler>(m, "MRMMapping", 
        R"doc(
A class to map targeted assays to chromatograms
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("mapExperiment", [](const OpenMS::MRMMapping& self, const OpenMS::MSExperiment& input_chromatograms, const OpenMS::TargetedExperiment& targeted_exp, OpenMS::MSExperiment& output) { return self.mapExperiment(input_chromatograms, targeted_exp, output); }, "input_chromatograms"_a, "targeted_exp"_a, "output"_a)
        .def("mapExperiment", [](const OpenMS::MRMMapping& self, const OpenMS::MSExperiment& input_chromatograms, const OpenSwath::LightTargetedExperiment& targeted_exp, OpenMS::MSExperiment& output) { return self.mapExperiment(input_chromatograms, targeted_exp, output); }, "input_chromatograms"_a, "targeted_exp"_a, "output"_a)
        ;

    // -----------------------------------------------------------------------
    // MRMTransitionGroupPicker
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMTransitionGroupPicker, OpenMS::DefaultParamHandler>(m, "MRMTransitionGroupPicker", 
        R"doc(
The MRMTransitionGroupPicker finds peaks in chromatograms that belong
to the same precursors
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("findLargestPeak", [](OpenMS::MRMTransitionGroupPicker& self, const std::vector<OpenMS::MSChromatogram>& picked_chroms) {
            int chr_idx = -1, peak_idx = -1;
            self.findLargestPeak(picked_chroms, chr_idx, peak_idx);
            return nb::make_tuple(chr_idx, peak_idx);
        }, "picked_chroms"_a, "Finds the largest peak across all chromatograms. Returns (chrom_idx, peak_idx)")
        .def("findWidestPeakIndices", [](const OpenMS::MRMTransitionGroupPicker& self, const std::vector<OpenMS::MSChromatogram>& picked_chroms) {
            OpenMS::Int chrom_idx = -1, point_idx = -1;
            self.findWidestPeakIndices(picked_chroms, chrom_idx, point_idx);
            return nb::make_tuple(chrom_idx, point_idx);
        }, "picked_chroms"_a, "Finds the widest peak across all chromatograms. Returns (chrom_idx, point_idx)")
        .def("pickTransitionGroup", [](OpenMS::MRMTransitionGroupPicker& self, OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& transition_group) { self.pickTransitionGroup(transition_group); }, "transition_group"_a, "Pick transition group")
        .def("remove_overlapping_features", [](OpenMS::MRMTransitionGroupPicker& self,
                std::vector<OpenMS::MSChromatogram>& picked_chroms,
                double best_left, double best_right) {
            self.remove_overlapping_features(picked_chroms, best_left, best_right);
        }, "picked_chroms"_a, "best_left"_a, "best_right"_a, "Remove overlapping features from picked chromatograms")
        .def("createMRMFeature", [](OpenMS::MRMTransitionGroupPicker& self,
                OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& transition_group,
                std::vector<OpenMS::MSChromatogram>& picked_chroms,
                const std::vector<OpenMS::MSChromatogram>& smoothed_chroms,
                int chr_idx, int peak_idx) {
            return self.createMRMFeature(transition_group, picked_chroms, smoothed_chroms, chr_idx, peak_idx);
        }, "transition_group"_a, "picked_chroms"_a, "smoothed_chroms"_a, "chr_idx"_a, "peak_idx"_a,
        "Create an MRMFeature from chromatograms and a specified peak")
        ;

    // -----------------------------------------------------------------------
    // MSDataAggregatingConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MSDataAggregatingConsumer, OpenMS::Interfaces::IMSDataConsumer>(m, "MSDataAggregatingConsumer", "Aggregates spectra by retention time")
        .def("__copy__", [](const OpenMS::MSDataAggregatingConsumer& self) { return OpenMS::MSDataAggregatingConsumer(self); })
        .def("__deepcopy__", [](const OpenMS::MSDataAggregatingConsumer& self, nb::dict) { return OpenMS::MSDataAggregatingConsumer(self); }, "memo"_a)
        .def("setExpectedSize", [](OpenMS::MSDataAggregatingConsumer& self, size_t p0, size_t p1) { return self.setExpectedSize(p0, p1); })
        .def("consumeSpectrum", [](OpenMS::MSDataAggregatingConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a)
        .def("consumeChromatogram", [](OpenMS::MSDataAggregatingConsumer& self, OpenMS::MSChromatogram& c) { return self.consumeChromatogram(c); }, "c"_a)
        .def("setExperimentalSettings", [](OpenMS::MSDataAggregatingConsumer& self, const OpenMS::ExperimentalSettings& p0) { return self.setExperimentalSettings(p0); })
        ;

    // -----------------------------------------------------------------------
    // MSDataSqlConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MSDataSqlConsumer, OpenMS::Interfaces::IMSDataConsumer>(m, "MSDataSqlConsumer", "A data consumer that inserts MS data into a SQLite database")
        .def(nb::init<OpenMS::String, size_t, int, bool, bool, double>())
        .def("flush", [](OpenMS::MSDataSqlConsumer& self) { return self.flush(); })
        .def("consumeSpectrum", [](OpenMS::MSDataSqlConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a, "Write a spectrum to the output file")
        .def("consumeChromatogram", [](OpenMS::MSDataSqlConsumer& self, OpenMS::MSChromatogram& c) { return self.consumeChromatogram(c); }, "c"_a, "Write a chromatogram to the output file")
        .def("setExpectedSize", [](OpenMS::MSDataSqlConsumer& self, size_t p0, size_t p1) { return self.setExpectedSize(p0, p1); })
        .def("setExperimentalSettings", [](OpenMS::MSDataSqlConsumer& self, const OpenMS::ExperimentalSettings& p0) { return self.setExperimentalSettings(p0); })
        ;

    // -----------------------------------------------------------------------
    // MSDataStoringConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MSDataStoringConsumer, OpenMS::Interfaces::IMSDataConsumer>(m, "MSDataStoringConsumer", 
        R"doc(
Consumer class that simply stores the data
This class is able to keep spectra and chromatograms passed to it in memory
and the data can be accessed through getData()
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MSDataStoringConsumer& self) { return OpenMS::MSDataStoringConsumer(self); })
        .def("__deepcopy__", [](const OpenMS::MSDataStoringConsumer& self, nb::dict) { return OpenMS::MSDataStoringConsumer(self); }, "memo"_a)
        .def("setExperimentalSettings", [](OpenMS::MSDataStoringConsumer& self, const OpenMS::ExperimentalSettings& settings) { return self.setExperimentalSettings(settings); }, "settings"_a, "Sets experimental settings")
        .def("setExpectedSize", [](OpenMS::MSDataStoringConsumer& self, size_t s_size, size_t c_size) { return self.setExpectedSize(s_size, c_size); }, "s_size"_a, "c_size"_a, "Sets expected size")
        .def("consumeSpectrum", [](OpenMS::MSDataStoringConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a)
        .def("consumeChromatogram", [](OpenMS::MSDataStoringConsumer& self, OpenMS::MSChromatogram& c) { return self.consumeChromatogram(c); }, "c"_a)
        .def("getData", [](const OpenMS::MSDataStoringConsumer& self) -> const OpenMS::MSExperiment & { return self.getData(); }, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // MSPFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MSPFile, OpenMS::DefaultParamHandler>(m, "MSPFile", "File adapter for MSP files (NIST spectra library)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MSPFile &>())
        .def("__copy__", [](const OpenMS::MSPFile& self) { return OpenMS::MSPFile(self); })
        .def("__deepcopy__", [](const OpenMS::MSPFile& self, nb::dict) { return OpenMS::MSPFile(self); }, "memo"_a)
        .def("load", [](OpenMS::MSPFile& self, const OpenMS::String& filename, OpenMS::PeptideIdentificationList& ids, OpenMS::MSExperiment& exp) { return self.load(filename, ids, exp); }, "filename"_a, "ids"_a, "exp"_a)
        .def("load", [](OpenMS::MSPFile& self, const OpenMS::String& filename, OpenMS::AnnotatedMSRun& annot_exp) { return self.load(filename, annot_exp); }, "filename"_a, "annot_exp"_a)
        .def("store", [](const OpenMS::MSPFile& self, const OpenMS::String& filename, const OpenMS::AnnotatedMSRun& exp) { return self.store(filename, exp); }, "filename"_a, "exp"_a, "Stores a map in a MSPFile file")
        ;

    // -----------------------------------------------------------------------
    // MassDecompositionAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MassDecompositionAlgorithm, OpenMS::DefaultParamHandler>(m, "MassDecompositionAlgorithm", 
        R"doc(
Mass decomposition algorithm, given a mass it suggests possible
compositions
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("getDecompositions", [](OpenMS::MassDecompositionAlgorithm& self, double weight) { std::vector<OpenMS::MassDecomposition> decomps; self.getDecompositions(decomps, weight); return decomps; }, "weight"_a, "Returns the possible decompositions given the weight")
        ;

    // -----------------------------------------------------------------------
    // MetaboliteFeatureDeconvolution
    // -----------------------------------------------------------------------
    auto metabolitefeaturedeconvolution_class = nb::class_<OpenMS::MetaboliteFeatureDeconvolution, OpenMS::DefaultParamHandler>(m, "MetaboliteFeatureDeconvolution", 
        R"doc(
DefaultParamHandler

An algorithm to decharge small molecule features (i.e. as found by FeatureFinder)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaboliteFeatureDeconvolution &>())
        .def("__copy__", [](const OpenMS::MetaboliteFeatureDeconvolution& self) { return OpenMS::MetaboliteFeatureDeconvolution(self); })
        .def("__deepcopy__", [](const OpenMS::MetaboliteFeatureDeconvolution& self, nb::dict) { return OpenMS::MetaboliteFeatureDeconvolution(self); }, "memo"_a)
        .def("compute", [](OpenMS::MetaboliteFeatureDeconvolution& self, const OpenMS::FeatureMap& fm_in, OpenMS::FeatureMap& fm_out, OpenMS::ConsensusMap& cons_map, OpenMS::ConsensusMap& cons_map_p) { return self.compute(fm_in, fm_out, cons_map, cons_map_p); }, "fm_in"_a, "fm_out"_a, "cons_map"_a, "cons_map_p"_a)
        ;
    nb::enum_<OpenMS::MetaboliteFeatureDeconvolution::CHARGEMODE>(metabolitefeaturedeconvolution_class, "CHARGEMODE",
        "Charge determination mode", nb::is_arithmetic())
        .value("QFROMFEATURE", OpenMS::MetaboliteFeatureDeconvolution::CHARGEMODE::QFROMFEATURE)
        .value("QHEURISTIC", OpenMS::MetaboliteFeatureDeconvolution::CHARGEMODE::QHEURISTIC)
        .value("QALL", OpenMS::MetaboliteFeatureDeconvolution::CHARGEMODE::QALL)
        ;

    // -----------------------------------------------------------------------
    // FeatureDeconvolution
    // -----------------------------------------------------------------------
    auto featuredeconvolution_class = nb::class_<OpenMS::FeatureDeconvolution, OpenMS::DefaultParamHandler>(m, "FeatureDeconvolution",
        "An algorithm to decharge features (i.e. as found by FeatureFinder)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureDeconvolution &>())
        .def("__copy__", [](const OpenMS::FeatureDeconvolution& self) { return OpenMS::FeatureDeconvolution(self); })
        .def("__deepcopy__", [](const OpenMS::FeatureDeconvolution& self, nb::dict) { return OpenMS::FeatureDeconvolution(self); }, "memo"_a)
        .def("compute", [](OpenMS::FeatureDeconvolution& self, const OpenMS::FeatureMap& fm_in, OpenMS::FeatureMap& fm_out, OpenMS::ConsensusMap& cons_map, OpenMS::ConsensusMap& cons_map_p) { self.compute(fm_in, fm_out, cons_map, cons_map_p); }, "fm_in"_a, "fm_out"_a, "cons_map"_a, "cons_map_p"_a)
        ;
    nb::enum_<OpenMS::FeatureDeconvolution::CHARGEMODE>(featuredeconvolution_class, "CHARGEMODE",
        "Charge determination mode", nb::is_arithmetic())
        .value("QFROMFEATURE", OpenMS::FeatureDeconvolution::CHARGEMODE::QFROMFEATURE)
        .value("QHEURISTIC", OpenMS::FeatureDeconvolution::CHARGEMODE::QHEURISTIC)
        .value("QALL", OpenMS::FeatureDeconvolution::CHARGEMODE::QALL)
        ;

    // -----------------------------------------------------------------------
    // MultiplexDeltaMassesGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MultiplexDeltaMassesGenerator, OpenMS::DefaultParamHandler>(m, "MultiplexDeltaMassesGenerator", 
        R"doc(
generates complete list of all possible mass shifts due to isotopic
labelling * * Isotopic labelling results in the shift of peptide
masses. * * For example in a Lys8/Arg10 SILAC labelled sample, some
peptides (the ones with one * Arg in their sequence) will show a
relative mass shift between light and heavy * partners of 10 Da. This
class constructs the complete list of all possible mass * shifts that
arise from isotopic labelling
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<OpenMS::String, int, std::map<OpenMS::String, double>>())
        .def("generateKnockoutDeltaMasses", [](OpenMS::MultiplexDeltaMassesGenerator& self) { return self.generateKnockoutDeltaMasses(); })
        .def("getDeltaMassesList", [](const OpenMS::MultiplexDeltaMassesGenerator& self) -> const std::vector<OpenMS::MultiplexDeltaMasses> & { return self.getDeltaMassesList(); }, nb::rv_policy::reference_internal)
        .def("getLabelShort", [](OpenMS::MultiplexDeltaMassesGenerator& self, const OpenMS::String& label) { return self.getLabelShort(label); }, "label"_a)
        .def("getLabelLong", [](OpenMS::MultiplexDeltaMassesGenerator& self, const OpenMS::String& label) { return self.getLabelLong(label); }, "label"_a)
        ;

    // -----------------------------------------------------------------------
    // NLargest
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::NLargest, OpenMS::DefaultParamHandler>(m, "NLargest", 
        R"doc(
DefaultParamHandler

NLargest removes all but the n largest peaks
)doc")
        .def(nb::init<>())
        .def(nb::init<unsigned int>())
        .def(nb::init<const OpenMS::NLargest &>())
        .def("__copy__", [](const OpenMS::NLargest& self) { return OpenMS::NLargest(self); })
        .def("__deepcopy__", [](const OpenMS::NLargest& self, nb::dict) { return OpenMS::NLargest(self); }, "memo"_a)
        .def("filterPeakSpectrum", [](OpenMS::NLargest& self, OpenMS::MSSpectrum& spectrum) { return self.filterPeakSpectrum(spectrum); }, "spectrum"_a, "Keep only n-largest peaks in spectrum")
        .def("filterPeakMap", [](OpenMS::NLargest& self, OpenMS::MSExperiment& exp) { return self.filterPeakMap(exp); }, "exp"_a, "Keep only n-largest peaks in each spectrum of a peak map")
        ;

    // -----------------------------------------------------------------------
    // Normalizer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Normalizer, OpenMS::DefaultParamHandler>(m, "Normalizer", 
        R"doc(
DefaultParamHandler

Normalizes the peak intensities spectrum-wise
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Normalizer &>())
        .def("__copy__", [](const OpenMS::Normalizer& self) { return OpenMS::Normalizer(self); })
        .def("__deepcopy__", [](const OpenMS::Normalizer& self, nb::dict) { return OpenMS::Normalizer(self); }, "memo"_a)
        .def("filterPeakSpectrum", [](const OpenMS::Normalizer& self, OpenMS::MSSpectrum& spectrum) { return self.filterPeakSpectrum(spectrum); }, "spectrum"_a, "Normalizes the peak spectrum")
        .def("filterPeakMap", [](const OpenMS::Normalizer& self, OpenMS::MSExperiment& exp) { return self.filterPeakMap(exp); }, "exp"_a, "Normalizes the peak map")
        ;

    // -----------------------------------------------------------------------
    // NucleicAcidSpectrumGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::NucleicAcidSpectrumGenerator, OpenMS::DefaultParamHandler>(m, "NucleicAcidSpectrumGenerator", 
        R"doc(
Generates theoretical spectra for nucleic acid sequences
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::NucleicAcidSpectrumGenerator &>())
        .def("__copy__", [](const OpenMS::NucleicAcidSpectrumGenerator& self) { return OpenMS::NucleicAcidSpectrumGenerator(self); })
        .def("__deepcopy__", [](const OpenMS::NucleicAcidSpectrumGenerator& self, nb::dict) { return OpenMS::NucleicAcidSpectrumGenerator(self); }, "memo"_a)
        .def("getSpectrum", [](const OpenMS::NucleicAcidSpectrumGenerator& self, const OpenMS::NASequence& oligo, int min_charge, int max_charge) { OpenMS::MSSpectrum spectrum; self.getSpectrum(spectrum, oligo, min_charge, max_charge); return spectrum; }, "oligo"_a, "min_charge"_a, "max_charge"_a, "Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters. If precursor_charge is set to 0 max_charge + 1 will be used")
        ;

    // -----------------------------------------------------------------------
    // OpenMSBuildInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Internal::OpenMSBuildInfo>(m, "OpenMSBuildInfo", "OpenMS class OpenMSBuildInfo")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Internal::OpenMSBuildInfo &>())
        .def("__copy__", [](const OpenMS::Internal::OpenMSBuildInfo& self) { return OpenMS::Internal::OpenMSBuildInfo(self); })
        .def("__deepcopy__", [](const OpenMS::Internal::OpenMSBuildInfo& self, nb::dict) { return OpenMS::Internal::OpenMSBuildInfo(self); }, "memo"_a)
        .def_static("isOpenMPEnabled", []() { return OpenMS::Internal::OpenMSBuildInfo::isOpenMPEnabled(); })
        .def_static("getBuildType", []() { return OpenMS::Internal::OpenMSBuildInfo::getBuildType(); })
        .def_static("getOpenMPMaxNumThreads", []() { return OpenMS::Internal::OpenMSBuildInfo::getOpenMPMaxNumThreads(); })
        .def_static("setOpenMPNumThreads", [](int num_threads) { return OpenMS::Internal::OpenMSBuildInfo::setOpenMPNumThreads(num_threads); }, "num_threads"_a)
        ;

    // -----------------------------------------------------------------------
    // OpenMSOSInfo
    // -----------------------------------------------------------------------
    auto openmsosinfo_class = nb::class_<OpenMS::Internal::OpenMSOSInfo>(m, "OpenMSOSInfo", "OpenMSOSInfo")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::Internal::OpenMSOSInfo& self) { return OpenMS::Internal::OpenMSOSInfo(self); })
        .def("__deepcopy__", [](const OpenMS::Internal::OpenMSOSInfo& self, nb::dict) { return OpenMS::Internal::OpenMSOSInfo(self); }, "memo"_a)
        .def("getOSAsString", [](const OpenMS::Internal::OpenMSOSInfo& self) { return self.getOSAsString(); })
        .def("getArchAsString", [](const OpenMS::Internal::OpenMSOSInfo& self) { return self.getArchAsString(); })
        .def("getOSVersionAsString", [](const OpenMS::Internal::OpenMSOSInfo& self) { return self.getOSVersionAsString(); })
        .def_static("getBinaryArchitecture", []() { return OpenMS::Internal::OpenMSOSInfo::getBinaryArchitecture(); })
        .def_static("getOSInfo", []() { return OpenMS::Internal::OpenMSOSInfo::getOSInfo(); })
        ;
    // OpenMS_OS enum nested under OpenMSOSInfo
    nb::enum_<OpenMS::Internal::OpenMS_OS>(openmsosinfo_class, "OpenMS_OS", nb::is_arithmetic())
        .value("OS_UNKNOWN", OpenMS::Internal::OpenMS_OS::OS_UNKNOWN)
        .value("OS_MACOS", OpenMS::Internal::OpenMS_OS::OS_MACOS)
        .value("OS_WINDOWS", OpenMS::Internal::OpenMS_OS::OS_WINDOWS)
        .value("OS_LINUX", OpenMS::Internal::OpenMS_OS::OS_LINUX)
        .value("SIZE_OF_OPENMS_OS", OpenMS::Internal::OpenMS_OS::SIZE_OF_OPENMS_OS)
        ;
    // OpenMS_Architecture enum nested under OpenMSOSInfo
    nb::enum_<OpenMS::Internal::OpenMS_Architecture>(openmsosinfo_class, "OpenMS_Architecture", nb::is_arithmetic())
        .value("ARCH_UNKNOWN", OpenMS::Internal::OpenMS_Architecture::ARCH_UNKNOWN)
        .value("ARCH_32BIT", OpenMS::Internal::OpenMS_Architecture::ARCH_32BIT)
        .value("ARCH_64BIT", OpenMS::Internal::OpenMS_Architecture::ARCH_64BIT)
        .value("SIZE_OF_OPENMS_ARCHITECTURE", OpenMS::Internal::OpenMS_Architecture::SIZE_OF_OPENMS_ARCHITECTURE)
        ;

    // -----------------------------------------------------------------------
    // PeakIntegrator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakIntegrator, OpenMS::DefaultParamHandler>(m, "PeakIntegrator", 
        R"doc(
DefaultParamHandler

Compute the area, background and shape metrics of a peak
The area computation is performed in integratePeak() and it supports
integration by simple sum of the intensity, integration by Simpson's rule
implementations for an odd number of unequally spaced points or integration
by the trapezoid rule
The background computation is performed in estimateBackground() and it
supports three different approaches to baseline correction, namely
computing a rectangular shape under the peak based on the minimum value of
the peak borders (vertical_division_min), a rectangular shape based on the
maximum value of the beak borders (vertical_division_max) or a trapezoidal
shape based on a straight line between the peak borders (base_to_base)
Peak shape metrics are computed in calculatePeakShapeMetrics() and multiple
metrics are supported
The containers supported by the methods are MSChromatogram and MSSpectrum
)doc")
        .def(nb::init<>())
        .def("integratePeak", [](const OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chromatogram, double left, double right) { return self.integratePeak(chromatogram, left, right); }, "chromatogram"_a, "left"_a, "right"_a)
        .def("integratePeak", [](const OpenMS::PeakIntegrator& self, const OpenMS::MSSpectrum& spectrum, double left, double right) { return self.integratePeak(spectrum, left, right); }, "spectrum"_a, "left"_a, "right"_a)
        .def("estimateBackground", [](const OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chromatogram, double left, double right, double peak_apex_pos) { return self.estimateBackground(chromatogram, left, right, peak_apex_pos); }, "chromatogram"_a, "left"_a, "right"_a, "peak_apex_pos"_a)
        .def("estimateBackground", [](const OpenMS::PeakIntegrator& self, const OpenMS::MSSpectrum& spectrum, double left, double right, double peak_apex_pos) { return self.estimateBackground(spectrum, left, right, peak_apex_pos); }, "spectrum"_a, "left"_a, "right"_a, "peak_apex_pos"_a)
        .def("calculatePeakShapeMetrics", [](const OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chromatogram, double left, double right, double peak_height, double peak_apex_pos) { return self.calculatePeakShapeMetrics(chromatogram, left, right, peak_height, peak_apex_pos); }, "chromatogram"_a, "left"_a, "right"_a, "peak_height"_a, "peak_apex_pos"_a)
        .def("calculatePeakShapeMetrics", [](const OpenMS::PeakIntegrator& self, const OpenMS::MSSpectrum& spectrum, double left, double right, double peak_height, double peak_apex_pos) { return self.calculatePeakShapeMetrics(spectrum, left, right, peak_height, peak_apex_pos); }, "spectrum"_a, "left"_a, "right"_a, "peak_height"_a, "peak_apex_pos"_a)

        .def("integratePeak", [](OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chrom, double left, double right) {
            return self.integratePeak(chrom, left, right);
        }, "chromatogram"_a, "left"_a, "right"_a, "Integrate peak in chromatogram")

        .def("integratePeak", [](OpenMS::PeakIntegrator& self, const OpenMS::MSSpectrum& spec, double left, double right) {
            return self.integratePeak(spec, left, right);
        }, "spectrum"_a, "left"_a, "right"_a, "Integrate peak in spectrum")

        .def("estimateBackground", [](OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chrom, double left, double right, double apex) {
            return self.estimateBackground(chrom, left, right, apex);
        }, "chromatogram"_a, "left"_a, "right"_a, "peak_apex_pos"_a, "Estimate background in chromatogram")

        .def("calculatePeakShapeMetrics", [](OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chrom, double left, double right, double height, double apex) {
            return self.calculatePeakShapeMetrics(chrom, left, right, height, apex);
        }, "chromatogram"_a, "left"_a, "right"_a, "peak_height"_a, "peak_apex_pos"_a, "Calculate peak shape metrics")
        .def("getDefaultParameters", [](OpenMS::PeakIntegrator& self) { OpenMS::Param params; self.getDefaultParameters(params); return params; }, "Returns default parameters")
        ;

    // -----------------------------------------------------------------------
    // PeakPickerChromatogram
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakPickerChromatogram, OpenMS::DefaultParamHandler>(m, "PeakPickerChromatogram", 
        R"doc(
DefaultParamHandler

The PeakPickerChromatogram finds peaks a single chromatogram
It uses the PeakPickerHiRes internally to find interesting seed candidates.
These candidates are then expanded and a right/left border of the peak is
searched
Additionally, overlapping peaks can be removed
)doc")
        .def(nb::init<>())
        .def("pickChromatogram", [](OpenMS::PeakPickerChromatogram& self, const OpenMS::MSChromatogram& chromatogram, OpenMS::MSChromatogram& picked_chrom) { return self.pickChromatogram(chromatogram, picked_chrom); }, "chromatogram"_a, "picked_chrom"_a)
        .def("pickChromatogram", [](OpenMS::PeakPickerChromatogram& self, const OpenMS::MSChromatogram& chromatogram, OpenMS::MSChromatogram& picked_chrom, OpenMS::MSChromatogram& smoothed_chrom) { return self.pickChromatogram(chromatogram, picked_chrom, smoothed_chrom); }, "chromatogram"_a, "picked_chrom"_a, "smoothed_chrom"_a)
        ;

    // -----------------------------------------------------------------------
    // PeakPickerIM
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakPickerIM, OpenMS::DefaultParamHandler>(m, "PeakPickerIM", "Peak picking algorithm for ion mobility data")
        .def(nb::init<>())
        .def("pickIMTraces", [](OpenMS::PeakPickerIM& self, OpenMS::MSSpectrum& spectrum) { return self.pickIMTraces(spectrum); }, "spectrum"_a, "Use trace detection for IM peak picking.")
        .def("pickIMCluster", [](const OpenMS::PeakPickerIM& self, OpenMS::MSSpectrum& spec) { return self.pickIMCluster(spec); }, "spec"_a, "Use clustering for IM peak picking.")
        .def("pickIMElutionProfiles", [](const OpenMS::PeakPickerIM& self, OpenMS::MSSpectrum& input) { return self.pickIMElutionProfiles(input); }, "input"_a, "Use elution profile detection for IM peak picking.")
        ;

    // -----------------------------------------------------------------------
    // PeptideAndProteinQuant
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideAndProteinQuant, OpenMS::DefaultParamHandler>(m, "PeptideAndProteinQuant", 
        R"doc(
Helper class for peptide and protein quantification based on feature
data annotated with IDs
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("readQuantData", [](OpenMS::PeptideAndProteinQuant& self, OpenMS::FeatureMap& features, const OpenMS::ExperimentalDesign& ed) { self.readQuantData(features, ed); }, "features"_a, "ed"_a,
            "Read quantitative data from a feature map")
        .def("readQuantData", [](OpenMS::PeptideAndProteinQuant& self, OpenMS::ConsensusMap& consensus, const OpenMS::ExperimentalDesign& ed) { self.readQuantData(consensus, ed); }, "consensus"_a, "ed"_a,
            "Read quantitative data from a consensus map")
        .def("quantifyPeptides", [](OpenMS::PeptideAndProteinQuant& self, const OpenMS::PeptideIdentificationList& peptides) { self.quantifyPeptides(peptides); },
            "peptides"_a = OpenMS::PeptideIdentificationList(),
            R"doc(
Compute peptide abundances.
Quantitative data must first be read via readQuantData().
Optional peptide-level protein inference information can be supplied via peptides.
)doc")
        .def("quantifyProteins", [](OpenMS::PeptideAndProteinQuant& self, const OpenMS::ProteinIdentification& proteins) { self.quantifyProteins(proteins); },
            "proteins"_a = OpenMS::ProteinIdentification(),
            R"doc(
Compute protein abundances.
Peptide abundances must be computed first with quantifyPeptides().
Optional protein inference information can be supplied via proteins.
)doc")
        .def("getStatistics", [](OpenMS::PeptideAndProteinQuant& self) -> const OpenMS::PeptideAndProteinQuant::Statistics & { return self.getStatistics(); }, nb::rv_policy::reference_internal,
            "Get summary statistics on quantification")
        .def("getPeptideResults", [](OpenMS::PeptideAndProteinQuant& self) { return self.getPeptideResults(); },
            "Get peptide abundance results as a dict mapping AASequence to PeptideData")
        .def("getProteinResults", [](OpenMS::PeptideAndProteinQuant& self) { return self.getProteinResults(); },
            "Get protein abundance results as a dict mapping protein accession (str) to ProteinData")
        ;

    // -----------------------------------------------------------------------
    // PosteriorErrorProbabilityModel
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::PosteriorErrorProbabilityModel, OpenMS::DefaultParamHandler>(m, "PosteriorErrorProbabilityModel", 
        R"doc(
Implements a mixture model of the inverse gumbel and the gauss
distribution or a gaussian mixture
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("fit", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, std::vector<double> search_engine_scores, const OpenMS::String& outlier_handling) { bool result = self.fit(search_engine_scores, outlier_handling); return nb::make_tuple(result, search_engine_scores); }, "search_engine_scores"_a, "outlier_handling"_a,
            R"doc(
Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables.
computeProbability can be used afterwards.
:param search_engine_scores: A vector which holds the data points
:return: Tuple of (success, sorted_scores). Note: scores are sorted from smallest to biggest value
)doc")
        .def("fit", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, std::vector<double> search_engine_scores, std::vector<double> probabilities, const OpenMS::String& outlier_handling) { bool result = self.fit(search_engine_scores, probabilities, outlier_handling); return nb::make_tuple(result, search_engine_scores, probabilities); }, "search_engine_scores"_a, "probabilities"_a, "outlier_handling"_a,
            R"doc(
Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
computeProbability can be used afterwards
Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
:param search_engine_scores: A vector which holds the data points
:param probabilities: Output vector for computed probabilities
:param outlier_handling: Outlier handling strategy
:return: Tuple of (success, sorted_scores, probabilities). success is `true` if algorithm has run through
)doc")
        .def("fillDensities", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, const std::vector<double>& x_scores) { std::vector<double> incorrect_density; std::vector<double> correct_density; self.fillDensities(x_scores, incorrect_density, correct_density); return nb::make_tuple(incorrect_density, correct_density); }, "x_scores"_a, "Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences")
        .def("fillLogDensities", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, const std::vector<double>& x_scores) { std::vector<double> incorrect_density; std::vector<double> correct_density; self.fillLogDensities(x_scores, incorrect_density, correct_density); return nb::make_tuple(incorrect_density, correct_density); }, "x_scores"_a, "Writes the log distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences")
        .def("computeLogLikelihood", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self, const std::vector<double>& incorrect_density, const std::vector<double>& correct_density) { return self.computeLogLikelihood(incorrect_density, correct_density); }, "incorrect_density"_a, "correct_density"_a, "Computes the Maximum Likelihood with a log-likelihood function")
        .def("pos_neg_mean_weighted_posteriors", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, const std::vector<double>& x_scores, const std::vector<double>& incorrect_posteriors) { return self.pos_neg_mean_weighted_posteriors(x_scores, incorrect_posteriors); }, "x_scores"_a, "incorrect_posteriors"_a)
        .def("getCorrectlyAssignedFitResult", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self) { return self.getCorrectlyAssignedFitResult(); }, "Returns estimated parameters for correctly assigned sequences. Fit should be used before")
        .def("getIncorrectlyAssignedFitResult", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self) { return self.getIncorrectlyAssignedFitResult(); }, "Returns estimated parameters for correctly assigned sequences. Fit should be used before")
        .def("getNegativePrior", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self) { return self.getNegativePrior(); }, "Returns the estimated negative prior probability")
        .def("computeProbability", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self, double score) { return self.computeProbability(score); }, "score"_a, "Returns the computed posterior error probability for a given score")
        .def("initPlots", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, std::vector<double> x_scores) { auto text_file = self.initPlots(x_scores); return nb::make_tuple(text_file, x_scores); }, "x_scores"_a, "Initializes the plots. Returns (TextFile, sorted_scores)")
        .def("getGumbelGnuplotFormula", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self, const OpenMS::Math::GaussFitter::GaussFitResult& params) { return self.getGumbelGnuplotFormula(params); }, "params"_a, "Returns the gnuplot formula of the fitted gumbel distribution")
        .def("getGaussGnuplotFormula", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self, const OpenMS::Math::GaussFitter::GaussFitResult& params) { return self.getGaussGnuplotFormula(params); }, "params"_a, "Returns the gnuplot formula of the fitted gauss distribution")
        .def("getBothGnuplotFormula", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self, const OpenMS::Math::GaussFitter::GaussFitResult& incorrect, const OpenMS::Math::GaussFitter::GaussFitResult& correct) { return self.getBothGnuplotFormula(incorrect, correct); }, "incorrect"_a, "correct"_a, "Returns the gnuplot formula of the fitted mixture distribution")
        .def("plotTargetDecoyEstimation", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, std::vector<double> target, std::vector<double> decoy) { self.plotTargetDecoyEstimation(target, decoy); }, "target"_a, "decoy"_a, "Plots the estimated distribution against target and decoy hits")
        .def("getSmallestScore", [](const OpenMS::Math::PosteriorErrorProbabilityModel& self) { return self.getSmallestScore(); }, "Returns the smallest score used in the last fit")
        .def("tryGnuplot", [](OpenMS::Math::PosteriorErrorProbabilityModel& self, const OpenMS::String& gp_file) { return self.tryGnuplot(gp_file); }, "gp_file"_a)
        ;

    // -----------------------------------------------------------------------
    // ProgressLogger
    // -----------------------------------------------------------------------
    auto progresslogger_class = nb::class_<OpenMS::ProgressLogger>(m, "ProgressLogger", 
        R"doc(
Base class for all classes that want to report their progress
Per default the progress log is disabled. Use setLogType to enable it
Use startProgress, setProgress and endProgress for the actual logging
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProgressLogger &>())
        .def("__copy__", [](const OpenMS::ProgressLogger& self) { return OpenMS::ProgressLogger(self); })
        .def("__deepcopy__", [](const OpenMS::ProgressLogger& self, nb::dict) { return OpenMS::ProgressLogger(self); }, "memo"_a)
        .def("setLogType", [](const OpenMS::ProgressLogger& self, OpenMS::ProgressLogger::LogType type) { return self.setLogType(type); }, "type"_a, "Sets the progress log that should be used. The default type is NONE!")
        .def("getLogType", [](const OpenMS::ProgressLogger& self) { return self.getLogType(); }, "Returns the type of progress log being used")
        .def("startProgress", [](const OpenMS::ProgressLogger& self, long begin, long end, const OpenMS::String& label) { return self.startProgress(begin, end, label); }, "begin"_a, "end"_a, "label"_a)
        .def("setProgress", [](const OpenMS::ProgressLogger& self, long value) { return self.setProgress(value); }, "value"_a, "Sets the current progress")
        .def("endProgress", [](const OpenMS::ProgressLogger& self, size_t bytes_processed) { return self.endProgress(bytes_processed); }, "bytes_processed"_a = 0, "Ends the progress display")
        .def("nextProgress", [](const OpenMS::ProgressLogger& self) { return self.nextProgress(); }, "Increment progress by 1 (according to range begin-end)")
        ;
    // LogType enum nested under ProgressLogger
    nb::enum_<OpenMS::ProgressLogger::LogType>(progresslogger_class, "LogType", "Enum for progress logging output type", nb::is_arithmetic())
        .value("CMD", OpenMS::ProgressLogger::LogType::CMD)
        .value("GUI", OpenMS::ProgressLogger::LogType::GUI)
        .value("NONE", OpenMS::ProgressLogger::LogType::NONE)

        .export_values();

    // -----------------------------------------------------------------------
    // AccurateMassSearchEngine
    // -----------------------------------------------------------------------
    auto accuratemasssearchengine_class = nb::class_<OpenMS::AccurateMassSearchEngine, OpenMS::DefaultParamHandler>(m, "AccurateMassSearchEngine", 
        R"doc(
An algorithm to search for exact mass matches from a spectrum against
a database (e.g. HMDB)
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("queryByMZ", [](const OpenMS::AccurateMassSearchEngine& self, const double& observed_mz, const int& observed_charge, const OpenMS::String& ion_mode, const OpenMS::EmpiricalFormula& observed_adduct) { std::vector<OpenMS::AccurateMassSearchResult> results; self.queryByMZ(observed_mz, observed_charge, ion_mode, results, observed_adduct); return results; }, "observed_mz"_a, "observed_charge"_a, "ion_mode"_a, "observed_adduct"_a)
        .def("queryByFeature", [](const OpenMS::AccurateMassSearchEngine& self, const OpenMS::Feature& feature, const size_t& feature_index, const OpenMS::String& ion_mode) { std::vector<OpenMS::AccurateMassSearchResult> results; self.queryByFeature(feature, feature_index, ion_mode, results); return results; }, "feature"_a, "feature_index"_a, "ion_mode"_a)
        .def("queryByConsensusFeature", [](const OpenMS::AccurateMassSearchEngine& self, const OpenMS::ConsensusFeature& cfeat, const size_t& cf_index, const size_t& number_of_maps, const OpenMS::String& ion_mode) { std::vector<OpenMS::AccurateMassSearchResult> results; self.queryByConsensusFeature(cfeat, cf_index, number_of_maps, ion_mode, results); return results; }, "cfeat"_a, "cf_index"_a, "number_of_maps"_a, "ion_mode"_a)
        .def("run", [](const OpenMS::AccurateMassSearchEngine& self, OpenMS::FeatureMap& p0, OpenMS::MzTab& p1) { return self.run(p0, p1); })
        .def("run", [](const OpenMS::AccurateMassSearchEngine& self, OpenMS::FeatureMap& p0, OpenMS::MzTabM& p1) { return self.run(p0, p1); })
        .def("run", [](const OpenMS::AccurateMassSearchEngine& self, OpenMS::ConsensusMap& p0, OpenMS::MzTab& p1) { return self.run(p0, p1); })
        .def("init", [](OpenMS::AccurateMassSearchEngine& self) { return self.init(); })
        ;
    def_ProgressLogger<OpenMS::AccurateMassSearchEngine>(accuratemasssearchengine_class);

    // -----------------------------------------------------------------------
    // AverageLinkage
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AverageLinkage, OpenMS::ProgressLogger>(m, "AverageLinkage", "AverageLinkage ClusterMethod")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::AverageLinkage &>())
        .def("__copy__", [](const OpenMS::AverageLinkage& self) { return OpenMS::AverageLinkage(self); })
        .def("__deepcopy__", [](const OpenMS::AverageLinkage& self, nb::dict) { return OpenMS::AverageLinkage(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // BaseGroupFinder
    // -----------------------------------------------------------------------
    auto basegroupfinder_class = nb::class_<OpenMS::BaseGroupFinder, OpenMS::DefaultParamHandler>(m, "BaseGroupFinder", 
        R"doc(
The base class of all element group finding algorithms
DefaultParamHandler
ProgressLogger
)doc")
        ;
    def_ProgressLogger<OpenMS::BaseGroupFinder>(basegroupfinder_class);

    // -----------------------------------------------------------------------
    // BasicProteinInferenceAlgorithm
    // -----------------------------------------------------------------------
    auto basicproteininferencealgorithm_class = nb::class_<OpenMS::BasicProteinInferenceAlgorithm, OpenMS::DefaultParamHandler>(m, "BasicProteinInferenceAlgorithm", 
        R"doc(
DefaultParamHandler
ProgressLogger

Algorithm class that implements simple protein inference by aggregation of peptide scores.
It has multiple parameter options like the aggregation method, when to distinguish peptidoforms,
and if you want to use shared peptides ("use_shared_peptides").
First, the best PSM per spectrum is used, then only the best PSM per peptidoform is aggregated.
Peptidoforms can optionally be distinguished via the treat_X_separate parameters:
- Modifications (modified sequence string)
- Charge states
The algorithm assumes posteriors or posterior error probabilities and converts to posteriors initially.
Possible aggregation methods that can be set via the parameter "aggregation_method" are:
- "best" (default)
- "sum"
- "product" (ignoring zeroes)
Annotation of the number of peptides used for aggregation can be disabled (see parameters).
Supports multiple runs but goes through them one by one iterating over the full PeptideIdentification vector.
Warning: Does not "link" the peptides to the resulting protein run. If you wish to do that you have to do
it manually.
Usage:
.. code-block:: python
from pyopenms import *
from urllib.request import urlretrieve
urlretrieve("https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/class_tests/openms/data/newMergerTest_out.idXML", "BasicProteinInference_test.idXML")
proteins = []
peptides = []
idf = IdXMLFile()
idf.load("BasicProteinInference_test.idXML", proteins, peptides);
bpia = BasicProteinInferenceAlgorithm()
p = bpia.getParameters();
p.setValue("min_peptides_per_protein", 0);
bpia.setParameters(p);
bpia.run(peptides, proteins);
#
hits = proteins[0].getHits()
print(hits[0].getScore()) # 0.6
print(hits[5].getScore()) # 0.9
print(hits[0].getMetaValue("nr_found_peptides")) # 1
print(hits[3].getMetaValue("nr_found_peptides")) # 2
)doc")
        .def(nb::init<>())
        .def("run", [](const OpenMS::BasicProteinInferenceAlgorithm& self, OpenMS::PeptideIdentificationList& pep_ids, std::vector<OpenMS::ProteinIdentification> prot_ids) {
            self.run(pep_ids, prot_ids);
            return prot_ids;
        }, "pep_ids"_a, "prot_ids"_a)
        .def("run", [](const OpenMS::BasicProteinInferenceAlgorithm& self, OpenMS::PeptideIdentificationList& pep_ids, OpenMS::ProteinIdentification& prot_id) { return self.run(pep_ids, prot_id); }, "pep_ids"_a, "prot_id"_a)
        .def("run", [](const OpenMS::BasicProteinInferenceAlgorithm& self, OpenMS::ConsensusMap& cmap, OpenMS::ProteinIdentification& prot_id, bool include_unassigned) { return self.run(cmap, prot_id, include_unassigned); }, "cmap"_a, "prot_id"_a, "include_unassigned"_a, 
            R"doc(
Performs basic aggregation-based inference on single ProteinIdentification run. See class help.
:param pep_ids: Vector of peptide identifications
:param prot_id: ProteinIdentification run with possible proteins. Scores will be overwritten and groups added.
:return: Writes its results into prot_ids
)doc")
        ;
    def_ProgressLogger<OpenMS::BasicProteinInferenceAlgorithm>(basicproteininferencealgorithm_class);
    // AggregationMethod enum nested under BasicProteinInferenceAlgorithm
    nb::enum_<OpenMS::BasicProteinInferenceAlgorithm::AggregationMethod>(basicproteininferencealgorithm_class, "AggregationMethod", nb::is_arithmetic())
        .value("PROD", OpenMS::BasicProteinInferenceAlgorithm::AggregationMethod::PROD)
        .value("SUM", OpenMS::BasicProteinInferenceAlgorithm::AggregationMethod::SUM)
        .value("BEST", OpenMS::BasicProteinInferenceAlgorithm::AggregationMethod::BEST)
        .export_values();

    // -----------------------------------------------------------------------
    // BayesianProteinInferenceAlgorithm
    // -----------------------------------------------------------------------
    auto bayesianproteininferencealgorithm_class = nb::class_<OpenMS::BayesianProteinInferenceAlgorithm, OpenMS::DefaultParamHandler>(m, "BayesianProteinInferenceAlgorithm", 
        R"doc(
DefaultParamHandler
ProgressLogger

Performs a Bayesian protein inference on Protein/Peptide identifications or ConsensusMap.
- Filters for best n PSMs per spectrum.
- Calculates and filters for best peptide per spectrum.
- Builds a k-partite graph from the structures.
- Finds and splits into connected components by DFS
- Extends the graph by adding layers from indist. protein groups, peptides with the same parents and optionally
some additional layers (peptide sequence, charge, replicate -> extended model = experimental)
- Builds a factor graph representation of a Bayesian network using the Evergreen library
See model param section. It is based on the Fido noisy-OR model with an option for
regularizing the number of proteins per peptide.
- Performs loopy belief propagation on the graph and queries protein, protein group and/or peptide posteriors
See loopy_belief_propagation param section.
- Learns best parameters via grid search if the parameters were not given in the param section.
- Writes posteriors to peptides and/or proteins and adds indistinguishable protein groups to the underlying
data structures.
- Can make use of OpenMP to parallelize over connected components.
Usage:
.. code-block:: python
from pyopenms import *
from urllib.request import urlretrieve
urlretrieve("https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/class_tests/openms/data/BayesianProteinInference_test.idXML", "BayesianProteinInference_test.idXML")
proteins = []
peptides = []
idf = IdXMLFile()
idf.load("BayesianProteinInference_test.idXML", proteins, peptides)
bpia = BayesianProteinInferenceAlgorithm()
p = bpia.getParameters()
p.setValue("update_PSM_probabilities", "false")
bpia.setParameters(p)
bpia.inferPosteriorProbabilities(proteins, peptides)
#
print(len(peptides)) # 9
print(peptides[0].getHits()[0].getScore()) # 0.6
print(proteins[0].getHits()[0].getScore()) # 0.624641
print(proteins[0].getHits()[1].getScore()) # 0.648346
)doc")
        .def(nb::init<unsigned int>())
        .def("inferPosteriorProbabilities", [](OpenMS::BayesianProteinInferenceAlgorithm& self, std::vector<OpenMS::ProteinIdentification> proteinIDs, OpenMS::PeptideIdentificationList& peptideIDs, bool greedy_group_resolution, std::optional<OpenMS::ExperimentalDesign> exp_des) {
            self.inferPosteriorProbabilities(proteinIDs, peptideIDs, greedy_group_resolution, exp_des);
            return proteinIDs;
        }, "proteinIDs"_a, "peptideIDs"_a, "greedy_group_resolution"_a, "exp_des"_a,
            R"doc(
Optionally adds indistinguishable protein groups with separate scores, too
Currently only takes first proteinID run and all peptides
:param proteinIDs: Vector of protein identifications
:param peptideIDs: Vector of peptide identifications
:returns: Updated proteinIDs with posterior probabilities
)doc")
        .def("inferPosteriorProbabilities", [](OpenMS::BayesianProteinInferenceAlgorithm& self, OpenMS::ConsensusMap& cmap, bool greedy_group_resolution, std::optional<OpenMS::ExperimentalDesign> exp_des) { return self.inferPosteriorProbabilities(cmap, greedy_group_resolution, exp_des); }, "cmap"_a, "greedy_group_resolution"_a, "exp_des"_a)
        ;
    def_ProgressLogger<OpenMS::BayesianProteinInferenceAlgorithm>(bayesianproteininferencealgorithm_class);

    // -----------------------------------------------------------------------
    // CachedMzMLHandler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Internal::CachedMzMLHandler, OpenMS::ProgressLogger>(m, "CachedMzMLHandler", 
        R"doc(
An class that uses on-disk caching to read and write spectra and
chromatograms
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("writeMemdump", [](const OpenMS::Internal::CachedMzMLHandler& self, const OpenMS::MSExperiment& exp, const OpenMS::String& out) { return self.writeMemdump(exp, out); }, "exp"_a, "out"_a, "Write complete spectra as a dump to the disk")
        .def("writeMetadata", [](OpenMS::Internal::CachedMzMLHandler& self, OpenMS::MSExperiment exp, const OpenMS::String& out_meta, bool addCacheMetaValue) { return self.writeMetadata(exp, out_meta, addCacheMetaValue); }, "exp"_a, "out_meta"_a, "addCacheMetaValue"_a = false, "Write only the meta data of an MSExperiment")
        .def("readMemdump", [](const OpenMS::Internal::CachedMzMLHandler& self, const OpenMS::String& filename) { OpenMS::MSExperiment exp_reading; self.readMemdump(exp_reading, filename); return exp_reading; }, "filename"_a, "Read all spectra from a dump from the disk")
        .def("createMemdumpIndex", [](OpenMS::Internal::CachedMzMLHandler& self, const OpenMS::String& filename) { return self.createMemdumpIndex(filename); }, "filename"_a, "Create an index on the location of all the spectra and chromatograms")
        .def("getSpectraIndex", [](const OpenMS::Internal::CachedMzMLHandler& self) -> const std::vector<std::streampos> & { return self.getSpectraIndex(); }, nb::rv_policy::reference_internal)
        .def("getChromatogramIndex", [](const OpenMS::Internal::CachedMzMLHandler& self) -> const std::vector<std::streampos> & { return self.getChromatogramIndex(); }, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // ChromatogramExtractor
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromatogramExtractor, OpenMS::ProgressLogger>(m, "ChromatogramExtractor", 
        R"doc(
The ChromatogramExtractor extracts chromatograms (intensity vs
retention time) from mass spectrometry data. * * This class provides
functionality to extract chromatographic traces from mass spectrometry
data * based on specified coordinates (m/z, retention time, and
optionally ion mobility values). * * The extractor supports two main
interfaces: * 1. Legacy interface: Takes a TargetedExperiment object
containing transitions and extracts * chromatograms at the m/z values
specified in those transitions. * 2. Modern interface: Takes a set of
ExtractionCoordinates that specify the exact coordinates * for
extraction. This provides more flexibility and control over the
extraction process. * The prepare_coordinates() helper function can
generate these coordinates for common * MS1 and MS2 extraction
scenarios. * * Key features: * - Supports both MS1 and MS2 level
extractions * - Configurable extraction window sizes in m/z dimension
(absolute or ppm) * - Multiple filter types available (Bartlett,
tophat) for signal processing * - Handles ion mobility data when
available * - Optimized for SWATH/DIA (Data Independent Acquisition)
experiments * - Progress logging capabilities through ProgressLogger
base class * * For MS2 extractions, the input data is expected to come
from a SWATH/DIA experiment * where precursor ions are fragmented in
wide isolation windows, allowing for extraction * of fragment ion
chromatograms. * * @see ChromatogramExtractorAlgorithm For the
underlying extraction algorithm implementation * @see
ExtractionCoordinates For the coordinate specification format
ProgressLogger
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChromatogramExtractor &>())
        .def("__copy__", [](const OpenMS::ChromatogramExtractor& self) { return OpenMS::ChromatogramExtractor(self); })
        .def("__deepcopy__", [](const OpenMS::ChromatogramExtractor& self, nb::dict) { return OpenMS::ChromatogramExtractor(self); }, "memo"_a)

        .def("extractChromatograms", [](OpenMS::ChromatogramExtractor& self,
                std::shared_ptr<OpenMS::SpectrumAccessOpenMS> input,
                nb::list output_py,
                std::vector<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates> extraction_coordinates,
                double mz_extraction_window,
                bool ppm,
                double im_extraction_window,
                const std::string& filter) {
            // Convert Python list to C++ vector
            std::vector<std::shared_ptr<OpenSwath::OSChromatogram>> output;
            for (size_t i = 0; i < nb::len(output_py); ++i) {
                output.push_back(nb::cast<std::shared_ptr<OpenSwath::OSChromatogram>>(output_py[i]));
            }
            self.extractChromatograms(input, output, extraction_coordinates, mz_extraction_window, ppm, im_extraction_window, filter);
            // Update the Python list in-place
            while (nb::len(output_py) > 0) { output_py.attr("pop")(); }
            for (auto& c : output) { output_py.append(nb::cast(c)); }
        }, "input"_a, "output"_a, "extraction_coordinates"_a, "mz_extraction_window"_a, "ppm"_a, "im_extraction_window"_a, "filter"_a)

        .def_static("prepare_coordinates", [](
                nb::list output_chromatograms_py,
                nb::list extraction_coordinates_py,
                OpenMS::TargetedExperiment& targeted,
                double rt_extraction_window,
                bool ms1,
                int ms1_isotopes) {
            std::vector<std::shared_ptr<OpenSwath::OSChromatogram>> output_chromatograms;
            std::vector<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates> extraction_coordinates;
            OpenMS::ChromatogramExtractor::prepare_coordinates(output_chromatograms, extraction_coordinates, targeted, rt_extraction_window, ms1, ms1_isotopes);
            // Update Python lists in-place
            while (nb::len(output_chromatograms_py) > 0) { output_chromatograms_py.attr("pop")(); }
            for (auto& c : output_chromatograms) { output_chromatograms_py.append(nb::cast(c)); }
            while (nb::len(extraction_coordinates_py) > 0) { extraction_coordinates_py.attr("pop")(); }
            for (auto& c : extraction_coordinates) { extraction_coordinates_py.append(nb::cast(c)); }
        }, "output_chromatograms"_a, "extraction_coordinates"_a, "targeted"_a, "rt_extraction_window"_a, "ms1"_a = false, "ms1_isotopes"_a = 0, "Prepare extraction coordinates from targeted experiment")
        ;

    // -----------------------------------------------------------------------
    // ChromatogramExtractorAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromatogramExtractorAlgorithm, OpenMS::ProgressLogger>(m, "ChromatogramExtractorAlgorithm", 
        R"doc(
The ChromatogramExtractorAlgorithm extracts chromatograms from a MS
data. * * It will take as input a set of transitions coordinates and
will extract * the signal of the provided map at the product ion m/z
and retention time * (rt) values specified by the extraction
coordinates. This interface only * expects a set of coordinates which
are up to the user to fill but a * convenient prepare_coordinates
function is provided (in the * ChromatogramExtractor class) to create
the coordinates for the most common * case of an MS2 and MS1
extraction. * * In the case of MS2 extraction, the map is assumed to
originate from a SWATH * (data-independent acquisition or DIA)
experiment. *
ProgressLogger
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChromatogramExtractorAlgorithm &>())
        .def("__copy__", [](const OpenMS::ChromatogramExtractorAlgorithm& self) { return OpenMS::ChromatogramExtractorAlgorithm(self); })
        .def("__deepcopy__", [](const OpenMS::ChromatogramExtractorAlgorithm& self, nb::dict) { return OpenMS::ChromatogramExtractorAlgorithm(self); }, "memo"_a)

        .def("extractChromatograms", [](OpenMS::ChromatogramExtractorAlgorithm& self,
                std::shared_ptr<OpenMS::SpectrumAccessOpenMS> input,
                std::vector<std::shared_ptr<OpenSwath::OSChromatogram>>& output,
                std::vector<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates> extraction_coordinates,
                double mz_extraction_window,
                bool ppm,
                double im_extraction_window,
                const std::string& filter) {
            self.extractChromatograms(input, output, extraction_coordinates, mz_extraction_window, ppm, im_extraction_window, filter);
        }, "input"_a, "output"_a, "extraction_coordinates"_a, "mz_extraction_window"_a, "ppm"_a, "im_extraction_window"_a, "filter"_a)
        ;

    // -----------------------------------------------------------------------
    // ConfidenceScoring
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConfidenceScoring, OpenMS::ProgressLogger>(m, "ConfidenceScoring", "OpenMS class ConfidenceScoring")
        .def(nb::init<bool>(), "test_mode"_a = false)
        .def("initialize", [](OpenMS::ConfidenceScoring& self, const OpenMS::TargetedExperiment& library, size_t n_decoys, size_t n_transitions, const OpenMS::TransformationDescription& rt_trafo) { return self.initialize(library, n_decoys, n_transitions, rt_trafo); }, "library"_a, "n_decoys"_a, "n_transitions"_a, "rt_trafo"_a)
        .def("initializeGlm", [](OpenMS::ConfidenceScoring& self, double intercept, double rt_coef, double int_coef) { return self.initializeGlm(intercept, rt_coef, int_coef); }, "intercept"_a, "rt_coef"_a, "int_coef"_a)
        .def("scoreMap", [](OpenMS::ConfidenceScoring& self, OpenMS::FeatureMap& features) { return self.scoreMap(features); }, "features"_a, "Score a feature map -> make sure the class is properly initialized")
        ;

    // -----------------------------------------------------------------------
    // DTA2DFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DTA2DFile, OpenMS::ProgressLogger>(m, "DTA2DFile", 
        R"doc(
DTA2D File adapter
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("load", [](OpenMS::DTA2DFile& self, const OpenMS::String& filename) {
            OpenMS::PeakMap exp;
            self.load(filename, exp);
            return exp;
        }, "filename"_a, "Loads a DTA2D file into an MSExperiment")
        .def("store", [](OpenMS::DTA2DFile& self, const OpenMS::String& filename, const OpenMS::PeakMap& map) { self.store(filename, map); }, "filename"_a, "map"_a, "Stores an MSExperiment to a DTA2D file")
        .def("getOptions", [](OpenMS::DTA2DFile& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal)
        .def("storeTIC", [](const OpenMS::DTA2DFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& map) { self.storeTIC(filename, map); }, "filename"_a, "map"_a, "Store TIC to file")
        ;

    // -----------------------------------------------------------------------
    // ElutionPeakDetection
    // -----------------------------------------------------------------------
    auto elutionpeakdetection_class = nb::class_<OpenMS::ElutionPeakDetection, OpenMS::DefaultParamHandler>(m, "ElutionPeakDetection", 
        R"doc(
Extracts chromatographic peaks from a mass trace
ProgressLogger
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("detectPeaks", [](OpenMS::ElutionPeakDetection& self, OpenMS::MassTrace& mt) { std::vector<OpenMS::MassTrace> single_mtraces; self.detectPeaks(mt, single_mtraces); return single_mtraces; }, "mt"_a)
        .def("detectPeaks", [](OpenMS::ElutionPeakDetection& self, std::vector<OpenMS::MassTrace> mt_vec) { std::vector<OpenMS::MassTrace> single_mtraces; self.detectPeaks(mt_vec, single_mtraces); return single_mtraces; }, "mt_vec"_a)
        .def("filterByPeakWidth", [](OpenMS::ElutionPeakDetection& self, std::vector<OpenMS::MassTrace> mt_vec) { std::vector<OpenMS::MassTrace> filtered; self.filterByPeakWidth(mt_vec, filtered); return filtered; }, "mt_vec"_a)
        .def("computeMassTraceNoise", [](OpenMS::ElutionPeakDetection& self, const OpenMS::MassTrace& p0) { return self.computeMassTraceNoise(p0); }, "Compute noise level (as RMSE of the actual signal and the smoothed signal)")
        .def("computeMassTraceSNR", [](OpenMS::ElutionPeakDetection& self, const OpenMS::MassTrace& p0) { return self.computeMassTraceSNR(p0); }, "Compute the signal to noise ratio (estimated by computeMassTraceNoise)")
        .def("computeApexSNR", [](OpenMS::ElutionPeakDetection& self, const OpenMS::MassTrace& p0) { return self.computeApexSNR(p0); }, "Compute the signal to noise ratio at the apex (estimated by computeMassTraceNoise)")
        .def("findLocalExtrema", [](const OpenMS::ElutionPeakDetection& self, const OpenMS::MassTrace& tr, const size_t& num_neighboring_peaks) { std::vector<size_t> chrom_maxes; std::vector<size_t> chrom_mins; self.findLocalExtrema(tr, num_neighboring_peaks, chrom_maxes, chrom_mins); return std::make_tuple(chrom_maxes, chrom_mins); }, "tr"_a, "num_neighboring_peaks"_a)
        .def("smoothData", [](const OpenMS::ElutionPeakDetection& self, OpenMS::MassTrace& mt, int win_size) { return self.smoothData(mt, win_size); }, "mt"_a, "win_size"_a)
        ;
    def_ProgressLogger<OpenMS::ElutionPeakDetection>(elutionpeakdetection_class);

    // -----------------------------------------------------------------------
    // FASTAFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FASTAFile, OpenMS::ProgressLogger>(m, "FASTAFile", 
        R"doc(
File adapter for FASTA files
Provides methods to load and store protein/peptide sequences in FASTA format.
Supports both batch loading (load/store) and streaming (readStart/readNext/writeStart/writeNext).
Usage:
.. code-block:: python
# Batch loading
entries = []
FASTAFile().load("proteins.fasta", entries)
for entry in entries:
print(entry.identifier, entry.sequence)
# Streaming (memory-efficient for large files)
ff = FASTAFile()
ff.readStart("proteins.fasta")
entry = FASTAEntry()
while ff.readNext(entry):
print(entry.identifier)
)doc")
        .def(nb::init<>())
        .def("readStart", [](OpenMS::FASTAFile& self, const OpenMS::String& filename) { return self.readStart(filename); }, "filename"_a)
        .def("atEnd", [](OpenMS::FASTAFile& self) { return self.atEnd(); }, "Boolean function to check if streams is at end of file")
        .def("writeStart", [](OpenMS::FASTAFile& self, const OpenMS::String& filename) { return self.writeStart(filename); }, "filename"_a)
        .def("writeEnd", [](OpenMS::FASTAFile& self) { return self.writeEnd(); }, "Closes the file (flush). Called implicitly when FASTAFile object does out of scope")
        .def("readNext", [](OpenMS::FASTAFile& self) {
            OpenMS::FASTAFile::FASTAEntry entry;
            bool success = self.readNext(entry);
            if (!success) return nb::make_tuple(nb::bool_(false), nb::cast(OpenMS::FASTAFile::FASTAEntry()));
            return nb::make_tuple(nb::bool_(true), nb::cast(entry));
        }, "Reads the next FASTA entry from file. Returns (success, FASTAEntry) tuple")
        .def("readNext", [](OpenMS::FASTAFile& self, OpenMS::FASTAFile::FASTAEntry& entry) {
            return self.readNext(entry);
        }, "entry"_a, "Reads the next FASTA entry from file into the given FASTAEntry. Returns True if entry was read, False if EOF")
        .def("writeNext", [](OpenMS::FASTAFile& self, const OpenMS::FASTAFile::FASTAEntry& entry) {
            self.writeNext(entry);
        }, "entry"_a, "Writes the given FASTAEntry to the file. Call writeStart() before and writeEnd() after")

        .def("load", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename, nb::list& output) {
            std::vector<OpenMS::FASTAFile::FASTAEntry> entries;
            {
                nb::gil_scoped_release release;
                self.load(filename, entries);
            }
            for (const auto& e : entries) {
                output.append(nb::cast(e));
            }
        }, "filename"_a, "data"_a, "Load a FASTA file. Appends FASTAEntry objects to output list")
        .def("load", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::FASTAFile::FASTAEntry> entries;
            {
                nb::gil_scoped_release release;
                self.load(filename, entries);
            }
            return entries;
        }, "filename"_a, "Load a FASTA file. Returns list of FASTAEntry objects")

        .def("store", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename, const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries) {
            nb::gil_scoped_release release;
            self.store(filename, entries);
        }, "filename"_a, "entries"_a, "Store a FASTA file. Takes list of FASTAEntry objects")
        .def("store", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename, const nb::list& entries) {
            std::vector<OpenMS::FASTAFile::FASTAEntry> fasta_entries;
            for (auto item : entries) {
                // Try to cast as FASTAEntry first, fall back to tuple
                try {
                    auto entry = nb::cast<OpenMS::FASTAFile::FASTAEntry>(item);
                    fasta_entries.push_back(entry);
                } catch (const nb::cast_error&) {
                    auto tup = nb::cast<nb::tuple>(item);
                    OpenMS::FASTAFile::FASTAEntry e;
                    e.identifier = nb::cast<std::string>(tup[0]);
                    e.description = nb::cast<std::string>(tup[1]);
                    e.sequence = nb::cast<std::string>(tup[2]);
                    fasta_entries.push_back(e);
                }
            }
            self.store(filename, fasta_entries);
        }, "filename"_a, "entries"_a, "Store a FASTA file. Takes list of FASTAEntry objects or (identifier, description, sequence) tuples")
        ;


    // -----------------------------------------------------------------------
    // FASTAEntry
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FASTAFile::FASTAEntry>(m, "FASTAEntry", "Represents a single FASTA entry with identifier, description, and sequence")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FASTAFile::FASTAEntry&>())
        .def("__copy__", [](const OpenMS::FASTAFile::FASTAEntry& self) { return OpenMS::FASTAFile::FASTAEntry(self); })
        .def("__deepcopy__", [](const OpenMS::FASTAFile::FASTAEntry& self, nb::dict) { return OpenMS::FASTAFile::FASTAEntry(self); }, "memo"_a)
        .def_rw("identifier", &OpenMS::FASTAFile::FASTAEntry::identifier, "Protein/sequence identifier")
        .def_rw("description", &OpenMS::FASTAFile::FASTAEntry::description, "Description from FASTA header")
        .def_rw("sequence", &OpenMS::FASTAFile::FASTAEntry::sequence, "The amino acid or nucleotide sequence")
        .def("headerMatches", &OpenMS::FASTAFile::FASTAEntry::headerMatches, "rhs"_a, "Returns True if identifier and description match the given entry")
        .def("sequenceMatches", &OpenMS::FASTAFile::FASTAEntry::sequenceMatches, "rhs"_a, "Returns True if the sequence matches the given entry")
        .def("__eq__", &OpenMS::FASTAFile::FASTAEntry::operator==)
        .def("__repr__", [](const OpenMS::FASTAFile::FASTAEntry& self) {
            std::ostringstream oss;
            oss << "FASTAEntry(identifier='" << self.identifier << "', description='" << self.description.substr(0, 50) << (self.description.size() > 50 ? "..." : "") << "')";
            return oss.str();
        })
        ;

    // -----------------------------------------------------------------------
    // FLASHDeconvAlgorithm
    // -----------------------------------------------------------------------
    auto flashdeconvalgorithm_class = nb::class_<OpenMS::FLASHDeconvAlgorithm, OpenMS::DefaultParamHandler>(m, "FLASHDeconvAlgorithm", 
        R"doc(
DefaultParamHandler
ProgressLogger

FLASHDeconv algorithm: ultrafast mass deconvolution algorithm for top down mass spectrometry dataset.
From MSSpectrum, this class outputs DeconvolvedSpectrum.
Deconvolution takes three steps:
1) decharging and select candidate masses - speed up via binning
2) collecting isotopes from the candidate masses and deisotoping - peak groups are defined here
3) scoring and filter out low scoring masses (i.e., peak groups)
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHDeconvAlgorithm &>())
        .def("__copy__", [](const OpenMS::FLASHDeconvAlgorithm& self) { return OpenMS::FLASHDeconvAlgorithm(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHDeconvAlgorithm& self, nb::dict) { return OpenMS::FLASHDeconvAlgorithm(self); }, "memo"_a)
        .def("getTolerances", [](const OpenMS::FLASHDeconvAlgorithm& self) { return self.getTolerances(); }, "Get calculated decoy averagine. Call after run() is called.")
        .def("run", [](OpenMS::FLASHDeconvAlgorithm& self, OpenMS::MSExperiment& map) { std::vector<OpenMS::DeconvolvedSpectrum> deconvolved_spectra; std::vector<OpenMS::FLASHHelperClasses::MassFeature> deconvolved_feature; { nb::gil_scoped_release release; self.run(map, deconvolved_spectra, deconvolved_feature); } return nb::make_tuple(deconvolved_spectra, deconvolved_feature); }, "map"_a)
        .def("getAveragine", [](OpenMS::FLASHDeconvAlgorithm& self) -> const OpenMS::FLASHHelperClasses::PrecalculatedAveragine & { return self.getAveragine(); }, nb::rv_policy::reference_internal, 
            R"doc(
Run FLASHDeconv algorithm for input_map and store deconvolved_spectra and deconvolved_features.
:param input_map: The input MSExperiment containing spectra to deconvolve
:param deconvolved_spectra: Output vector to store deconvolved spectra
:param deconvolved_features: Output vector to store mass features
Averagine access
)doc")
        .def("getDecoyAveragine", [](OpenMS::FLASHDeconvAlgorithm& self) -> const OpenMS::FLASHHelperClasses::PrecalculatedAveragine & { return self.getDecoyAveragine(); }, nb::rv_policy::reference_internal, "Get calculated averagine. Call after run() is called.")
        .def("getNoiseDecoyWeight", [](const OpenMS::FLASHDeconvAlgorithm& self) { return self.getNoiseDecoyWeight(); }, "Get mass tolerances per MS level.")
        .def_static("getScanNumber", [](const OpenMS::MSExperiment& map, size_t index) { return OpenMS::FLASHDeconvAlgorithm::getScanNumber(map, index); }, "map"_a, "index"_a, "Get noise decoy weight determined during q-value calculation.")
        ;
    def_ProgressLogger<OpenMS::FLASHDeconvAlgorithm>(flashdeconvalgorithm_class);

    // -----------------------------------------------------------------------
    // FeatureFinderAlgorithmPicked
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureFinderAlgorithmPicked, OpenMS::DefaultParamHandler>(m, "FeatureFinderAlgorithmPicked", 
        R"doc(
FeatureFinder algorithm for picked (centroided) data
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::FeatureFinderAlgorithmPicked& self, OpenMS::MSExperiment& input_map, OpenMS::FeatureMap& features, const OpenMS::Param& param, const OpenMS::FeatureMap& seeds) { nb::gil_scoped_release release; return self.run(input_map, features, param, seeds); }, "input_map"_a, "features"_a, "param"_a, "seeds"_a)
        ;

    // -----------------------------------------------------------------------
    // FeatureFinderMultiplexAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureFinderMultiplexAlgorithm, OpenMS::DefaultParamHandler>(m, "FeatureFinderMultiplexAlgorithm", 
        R"doc(
FeatureFinderMultiplexAlgorithm is a tool for the fully automated
analysis of quantitative proteomics data. It detects pairs of isotopic
envelopes with fixed m/z separation. It requires no prior sequence
identification of the peptides and works on both profile or centroided
spectra. In what follows we outline the algorithm
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::FeatureFinderMultiplexAlgorithm& self, OpenMS::MSExperiment& exp, bool progress) { nb::gil_scoped_release release; return self.run(exp, progress); }, "exp"_a, "progress"_a, "Main method for feature detection")
        .def("getFeatureMap", [](OpenMS::FeatureFinderMultiplexAlgorithm& self) -> OpenMS::FeatureMap & { return self.getFeatureMap(); }, nb::rv_policy::reference_internal)
        .def("getConsensusMap", [](OpenMS::FeatureFinderMultiplexAlgorithm& self) -> OpenMS::ConsensusMap & { return self.getConsensusMap(); }, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // FeatureFindingMetabo
    // -----------------------------------------------------------------------
    auto featurefindingmetabo_class = nb::class_<OpenMS::FeatureFindingMetabo, OpenMS::DefaultParamHandler>(m, "FeatureFindingMetabo", 
        R"doc(
Internal structure used in @ref FeatureFindingMetabo that keeps track
of a feature hypothesis (isotope group hypothesis)
ProgressLogger
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::FeatureFindingMetabo& self, std::vector<OpenMS::MassTrace> input_mtraces, OpenMS::FeatureMap& output_featmap) { std::vector<std::vector<OpenMS::MSChromatogram>> output_chromatograms; { nb::gil_scoped_release release; self.run(input_mtraces, output_featmap, output_chromatograms); } return nb::make_tuple(input_mtraces, output_chromatograms); }, "input_mtraces"_a, "output_featmap"_a)
        ;
    def_ProgressLogger<OpenMS::FeatureFindingMetabo>(featurefindingmetabo_class);

    // -----------------------------------------------------------------------
    // FeatureGroupingAlgorithmKD
    // -----------------------------------------------------------------------
    auto featuregroupingalgorithmkd_class = nb::class_<OpenMS::FeatureGroupingAlgorithmKD, OpenMS::FeatureGroupingAlgorithm>(m, "FeatureGroupingAlgorithmKD", 
        R"doc(
Proxy for a (potential) cluster
FeatureGroupingAlgorithm
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("group", [](OpenMS::FeatureGroupingAlgorithmKD& self, const std::vector<OpenMS::FeatureMap>& maps, OpenMS::ConsensusMap& out) { self.group(maps, out); }, "maps"_a, "out"_a)
        .def("group", [](OpenMS::FeatureGroupingAlgorithmKD& self, const std::vector<OpenMS::ConsensusMap>& maps, OpenMS::ConsensusMap& out) { self.group(maps, out); }, "maps"_a, "out"_a, "Group consensus maps")
        ;
    def_ProgressLogger<OpenMS::FeatureGroupingAlgorithmKD>(featuregroupingalgorithmkd_class);

    // -----------------------------------------------------------------------
    // GNPSMGFFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::GNPSMGFFile, OpenMS::DefaultParamHandler>(m, "GNPSMGFFile", "OpenMS class GNPSMGFFile")
        .def(nb::init<>())
        .def("store", [](const OpenMS::GNPSMGFFile& self, const OpenMS::String& consensus_file_path, const std::vector<OpenMS::String>& mzml_file_paths, const OpenMS::String& out) { return self.store(consensus_file_path, mzml_file_paths, out); }, "consensus_file_path"_a, "mzml_file_paths"_a, "out"_a, "Export consensus file from default workflow to GNPS MGF format")
        ;

    // -----------------------------------------------------------------------
    // GaussFilter
    // -----------------------------------------------------------------------
    auto gaussfilter_class = nb::class_<OpenMS::GaussFilter, OpenMS::ProgressLogger>(m, "GaussFilter", "OpenMS class GaussFilter")
        .def(nb::init<>())
        .def("filter", [](OpenMS::GaussFilter& self, OpenMS::MSSpectrum& spectrum) { return self.filter(spectrum); }, "spectrum"_a, "Smoothes an MSSpectrum containing profile data")
        .def("filter", [](OpenMS::GaussFilter& self, OpenMS::MSChromatogram& chromatogram) { return self.filter(chromatogram); }, "chromatogram"_a, "Smoothes an MSSpectrum containing profile data")
        .def("filter", [](OpenMS::GaussFilter& self, OpenMS::Mobilogram& mobilogram) { return self.filter(mobilogram); }, "mobilogram"_a, "Smoothes an MSSpectrum containing profile data")
        .def("filterExperiment", [](OpenMS::GaussFilter& self, OpenMS::MSExperiment& map) { return self.filterExperiment(map); }, "map"_a, "Smoothes an MSExperiment containing profile data")
        ;
    def_DefaultParamHandler<OpenMS::GaussFilter>(gaussfilter_class);

    // -----------------------------------------------------------------------
    // ModifiedSincSmoother
    // -----------------------------------------------------------------------
    auto modifiedsincsmoother_class = nb::class_<OpenMS::ModifiedSincSmoother, OpenMS::ProgressLogger>(m, "ModifiedSincSmoother",
        R"doc(
Modified sinc smoother for profile data.
Two variants: MS (better stopband, larger kernel) and MS1 (smaller kernel).
Based on Schmid, Rath & Diebold, ACS Meas. Sci. Au 2022.
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def(nb::init<bool, int, int>(), "isMS1"_a, "degree"_a, "m"_a,
            "Create smoother. isMS1: true for MS1 variant (smaller kernel), false for MS variant.")
        .def("filter", [](OpenMS::ModifiedSincSmoother& self, OpenMS::MSSpectrum& spectrum) {
            return self.filter(spectrum); }, "spectrum"_a, "Smooth an MSSpectrum in-place")
        .def("filter", [](OpenMS::ModifiedSincSmoother& self, OpenMS::MSChromatogram& chromatogram) {
            return self.filter(chromatogram); }, "chromatogram"_a, "Smooth an MSChromatogram in-place")
        .def("filter", [](OpenMS::ModifiedSincSmoother& self, OpenMS::Mobilogram& mobilogram) {
            return self.filter(mobilogram); }, "mobilogram"_a, "Smooth a Mobilogram in-place")
        .def("filterExperiment", [](OpenMS::ModifiedSincSmoother& self, OpenMS::MSExperiment& map) {
            return self.filterExperiment(map); }, "map"_a, "Smooth all spectra and chromatograms in a PeakMap")
        .def("smooth", [](OpenMS::ModifiedSincSmoother& self, const std::vector<double>& data) {
            return self.smooth(data); }, "data"_a, "Smooth a vector of intensity values, returns smoothed vector")
        .def_static("bandwidthToM", &OpenMS::ModifiedSincSmoother::bandwidthToM,
            "isMS1"_a, "degree"_a, "bandwidth"_a, "Convert frequency bandwidth to kernel half-width m")
        .def_static("noiseGainToM", &OpenMS::ModifiedSincSmoother::noiseGainToM,
            "isMS1"_a, "degree"_a, "noiseGain"_a, "Convert noise gain to kernel half-width m")
        .def_static("savitzkyGolayBandwidth", &OpenMS::ModifiedSincSmoother::savitzkyGolayBandwidth,
            "degree"_a, "m"_a, "Compute equivalent Savitzky-Golay bandwidth")
        ;
    def_DefaultParamHandler<OpenMS::ModifiedSincSmoother>(modifiedsincsmoother_class);

    // -----------------------------------------------------------------------
    // InternalCalibration
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::InternalCalibration, OpenMS::ProgressLogger>(m, "InternalCalibration", 
        R"doc(
A mass recalibration method using linear/quadratic interpolation
(robust/weighted) of given reference masses
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("fillCalibrants", [](OpenMS::InternalCalibration& self, const OpenMS::PeakMap& exp,
                                   const std::vector<OpenMS::InternalCalibration::LockMass>& ref_masses,
                                   double tol_ppm, bool lock_require_mono, bool lock_require_iso, bool verbose) {
            OpenMS::CalibrationData failed_lock_masses;
            const OpenMS::Size found = self.fillCalibrants(exp, ref_masses, tol_ppm, lock_require_mono, lock_require_iso, failed_lock_masses, verbose);
            return nb::make_tuple(found, failed_lock_masses);
        }, "exp"_a, "ref_masses"_a, "tol_ppm"_a, "lock_require_mono"_a, "lock_require_iso"_a, "verbose"_a = true,
           "Extract calibrants from raw spectra. Returns (number_found, failed_lock_masses)")
        .def("fillCalibrants", [](OpenMS::InternalCalibration& self, const OpenMS::FeatureMap& fm, double tol_ppm) {
            return self.fillCalibrants(fm, tol_ppm);
        }, "fm"_a, "tol_ppm"_a, "Extract calibrants from feature identifications")
        .def("fillCalibrants", [](OpenMS::InternalCalibration& self, const OpenMS::PeptideIdentificationList& pep_ids, double tol_ppm) {
            return self.fillCalibrants(pep_ids, tol_ppm);
        }, "pep_ids"_a, "tol_ppm"_a, "Extract calibrants from peptide identifications")
        .def("getCalibrationPoints", [](const OpenMS::InternalCalibration& self) -> const OpenMS::CalibrationData & { return self.getCalibrationPoints(); }, nb::rv_policy::reference_internal, 
            R"doc(
Extract calibrants from identifications\n
Extracts only the first hit from each peptide identification
Hits are sorted beforehand
Ambiguities should be resolved before, e.g. using IDFilter\n
Unassigned peptide identifications are also taken into account!
RT and m/z are naturally taken from the IDs, since to feature is assigned
If you do not want these IDs, remove them from the feature map before calling this function\n
A filtering step is done in the m/z dimension using 'tol_ppm'
Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
larger outliers are removed before accepting an ID as calibrant
:param pep_ids: Peptide ids (e.g. from an idXML file)
:param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
:return: Number of calibration masses found
)doc")
        .def_static("applyTransformation", [](std::vector<OpenMS::Precursor> pcs, const OpenMS::MZTrafoModel& trafo) {
            OpenMS::InternalCalibration::applyTransformation(pcs, trafo);
            return pcs;
        }, "pcs"_a, "trafo"_a,
            R"doc(
Apply calibration to data\n
For each spectrum, a calibration model will be computed and applied.
Make sure to call fillCalibrants() before, so a model can be created.\n
The MSExperiment will be sorted by RT and m/z if unsorted.
:param exp: MSExperiment holding the Raw data to calibrate
:param target_mslvl: MS-levels where calibration should be applied to
:param model_type: Linear or quadratic model; select based on your instrument
:param rt_chunk: RT-window size (one-sided) of calibration points to collect around each spectrum. Set to negative values, to build one global model instead.
:param use_RANSAC: Remove outliers before fitting a model?!
:param post_ppm_median: The median ppm error of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
:param post_ppm_MAD: The median absolute deviation of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
:param file_models: Output CSV filename, where model parameters are written to (pass empty string to skip)
:param file_models_plot: Output PNG image model parameters (pass empty string to skip)
:param file_residuals: Output CSV filename, where ppm errors of calibrants before and after model fitting parameters are written to (pass empty string to skip)
:param file_residuals_plot: Output PNG image of the ppm errors of calibrants (pass empty string to skip)
:param rscript_executable: Full path to the Rscript executable
:returns: Updated precursors with calibration applied
)doc")
        .def("calibrate", [](OpenMS::InternalCalibration& self, OpenMS::PeakMap& exp, const OpenMS::IntList& target_mslvl,
                              OpenMS::MZTrafoModel::MODELTYPE model_type, double rt_chunk, bool use_RANSAC,
                              double post_ppm_median, double post_ppm_MAD,
                              const OpenMS::String& file_models, const OpenMS::String& file_models_plot,
                              const OpenMS::String& file_residuals, const OpenMS::String& file_residuals_plot,
                              const OpenMS::String& rscript_executable) {
            return self.calibrate(exp, target_mslvl, model_type, rt_chunk, use_RANSAC,
                                  post_ppm_median, post_ppm_MAD, file_models, file_models_plot,
                                  file_residuals, file_residuals_plot, rscript_executable);
        }, "exp"_a, "target_mslvl"_a, "model_type"_a, "rt_chunk"_a, "use_RANSAC"_a,
           "post_ppm_median"_a, "post_ppm_MAD"_a, "file_models"_a = "", "file_models_plot"_a = "",
           "file_residuals"_a = "", "file_residuals_plot"_a = "", "rscript_executable"_a = "Rscript",
           "Calibrate an experiment using previously collected calibration points")
        .def_static("applyTransformation", [](OpenMS::MSSpectrum& spec, const std::vector<int>& target_mslvl, const OpenMS::MZTrafoModel& trafo) { return OpenMS::InternalCalibration::applyTransformation(spec, target_mslvl, trafo); }, "spec"_a, "target_mslvl"_a, "trafo"_a)
        .def_static("applyTransformation", [](OpenMS::MSExperiment& exp, const std::vector<int>& target_mslvl, const OpenMS::MZTrafoModel& trafo) { return OpenMS::InternalCalibration::applyTransformation(exp, target_mslvl, trafo); }, "exp"_a, "target_mslvl"_a, "trafo"_a)
        ;

    // -----------------------------------------------------------------------
    // LabeledPairFinder
    // -----------------------------------------------------------------------
    auto labeledpairfinder_class = nb::class_<OpenMS::LabeledPairFinder, OpenMS::BaseGroupFinder>(m, "LabeledPairFinder", 
        R"doc(
The LabeledPairFinder allows the matching of labeled features (features with a fixed distance)
Finds feature pairs that have a defined distance in RT and m/z in the same map
BaseGroupFinder
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::LabeledPairFinder& self, const std::vector<OpenMS::ConsensusMap>& input_maps, OpenMS::ConsensusMap& result_map) { self.run(input_maps, result_map); }, "input_maps"_a, "result_map"_a)
        ;
    def_ProgressLogger<OpenMS::LabeledPairFinder>(labeledpairfinder_class);

    // -----------------------------------------------------------------------
    // LinearResamplerAlign
    // -----------------------------------------------------------------------
    auto linearresampleralign_class = nb::class_<OpenMS::LinearResamplerAlign, OpenMS::DefaultParamHandler>(m, "LinearResamplerAlign",
        R"doc(
Linear Resampling of raw data with alignment
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::LinearResamplerAlign &>())
        .def("__copy__", [](const OpenMS::LinearResamplerAlign& self) { return OpenMS::LinearResamplerAlign(self); })
        .def("__deepcopy__", [](const OpenMS::LinearResamplerAlign& self, nb::dict) { return OpenMS::LinearResamplerAlign(self); }, "memo"_a)
        .def("raster", [](OpenMS::LinearResamplerAlign& self, OpenMS::MSSpectrum& spectrum) { return self.raster(spectrum); }, "spectrum"_a, "Applies the resampling algorithm to an MSSpectrum")
        .def("raster", [](OpenMS::LinearResamplerAlign& self, OpenMS::MSChromatogram& chromatogram) { return self.raster(chromatogram); }, "chromatogram"_a, "Applies the resampling algorithm to an MSChromatogram")
        .def("raster_align", [](OpenMS::LinearResamplerAlign& self, OpenMS::MSSpectrum& spectrum, double start_pos, double end_pos) { return self.raster_align(spectrum, start_pos, end_pos); }, "spectrum"_a, "start_pos"_a, "end_pos"_a, "Resamples an MSSpectrum onto an explicit raster [start_pos, end_pos]")
        .def("raster_align", [](OpenMS::LinearResamplerAlign& self, OpenMS::MSChromatogram& chromatogram, double start_pos, double end_pos) { return self.raster_align(chromatogram, start_pos, end_pos); }, "chromatogram"_a, "start_pos"_a, "end_pos"_a, "Resamples an MSChromatogram onto an explicit raster [start_pos, end_pos]")
        .def("rasterExperiment", [](OpenMS::LinearResamplerAlign& self, OpenMS::MSExperiment& exp) { return self.rasterExperiment(exp); }, "exp"_a, "Resamples all spectra in an MSExperiment")
        ;
    def_ProgressLogger<OpenMS::LinearResamplerAlign>(linearresampleralign_class);

    // -----------------------------------------------------------------------
    // MRMAssay
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMAssay, OpenMS::ProgressLogger>(m, "MRMAssay", 
        R"doc(
Generate assays from a TargetedExperiment
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("restrictTransitions", [](OpenMS::MRMAssay& self, OpenMS::TargetedExperiment& exp, double lower_mz_limit, double upper_mz_limit, const std::vector<std::pair<double, double>>& swathes) { return self.restrictTransitions(exp, lower_mz_limit, upper_mz_limit, swathes); }, "exp"_a, "lower_mz_limit"_a, "upper_mz_limit"_a, "swathes"_a, 
            R"doc(
Annotates and filters transitions in a TargetedExperiment
:param exp: The input, unfiltered transitions
:param precursor_mz_threshold: The precursor m/z threshold in Th for annotation
:param product_mz_threshold: The product m/z threshold in Th for annotation
:param fragment_types: The fragment types to consider for annotation
:param fragment_charges: The fragment charges to consider for annotation
:param enable_specific_losses: Whether specific neutral losses should be considered
:param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
:param round_decPow: Round product m/z values to decimal power (default: -4)
)doc")
        .def("detectingTransitions", [](OpenMS::MRMAssay& self, OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions) { return self.detectingTransitions(exp, min_transitions, max_transitions); }, "exp"_a, "min_transitions"_a, "max_transitions"_a, 
            R"doc(
Restrict and filter transitions in a TargetedExperiment
:param exp: The input, unfiltered transitions
:param lower_mz_limit: The lower product m/z limit in Th
:param upper_mz_limit: The upper product m/z limit in Th
:param swathes: The swath window settings (to exclude fragment ions falling into the precursor isolation window)
)doc")
        .def("filterMinMaxTransitionsCompound", [](OpenMS::MRMAssay& self, OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions) { return self.filterMinMaxTransitionsCompound(exp, min_transitions, max_transitions); }, "exp"_a, "min_transitions"_a, "max_transitions"_a, 
            R"doc(
Select detecting fragment ions
:param exp: The input, unfiltered transitions
:param min_transitions: The minimum number of transitions required per assay
:param max_transitions: The maximum number of transitions required per assay
)doc")
        .def("filterUnreferencedDecoysCompound", [](OpenMS::MRMAssay& self, OpenMS::TargetedExperiment& exp) { return self.filterUnreferencedDecoysCompound(exp); }, "exp"_a, 
            R"doc(
Filters target and decoy transitions by intensity, only keeping the top N transitions
:param exp: The transition list which will be filtered
:param min_transitions: The minimum number of transitions required per assay (targets only)
:param max_transitions: The maximum number of transitions allowed per assay
)doc")
        .def("reannotateTransitions", [](OpenMS::MRMAssay& self, OpenMS::TargetedExperiment& exp, double precursor_mz_threshold, double product_mz_threshold, const std::vector<OpenMS::String>& fragment_types, const std::vector<size_t>& fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow) { self.reannotateTransitions(exp, precursor_mz_threshold, product_mz_threshold, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, round_decPow); }, "exp"_a, "precursor_mz_threshold"_a, "product_mz_threshold"_a, "fragment_types"_a, "fragment_charges"_a, "enable_specific_losses"_a, "enable_unspecific_losses"_a, "round_decPow"_a = -4)
        .def("uisTransitions", [](OpenMS::MRMAssay& self, OpenMS::TargetedExperiment& exp, const std::vector<OpenMS::String>& fragment_types, const std::vector<size_t>& fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, bool enable_ms2_precursors, double mz_threshold, const std::vector<std::pair<double, double>>& swathes, int round_decPow, size_t max_num_alternative_localizations, int shuffle_seed, bool disable_decoy_transitions) { self.uisTransitions(exp, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, enable_ms2_precursors, mz_threshold, swathes, round_decPow, max_num_alternative_localizations, shuffle_seed, disable_decoy_transitions); }, "exp"_a, "fragment_types"_a, "fragment_charges"_a, "enable_specific_losses"_a, "enable_unspecific_losses"_a, "enable_ms2_precursors"_a, "mz_threshold"_a, "swathes"_a, "round_decPow"_a = -4, "max_num_alternative_localizations"_a = 20, "shuffle_seed"_a = -1, "disable_decoy_transitions"_a = false)
        ;

    // -----------------------------------------------------------------------
    // MRMDecoy
    // -----------------------------------------------------------------------
    auto mrmdecoy_class = nb::class_<OpenMS::MRMDecoy, OpenMS::DefaultParamHandler>(m, "MRMDecoy", 
        R"doc(
This class generates a TargetedExperiment object with decoys based on
a TargetedExperiment object
ProgressLogger
)doc")
        .def(nb::init<>())
        .def_static("findFixedResidues", [](const OpenMS::String& sequence, bool keepN, bool keepC, const OpenMS::String& keep_const_pattern) { return OpenMS::MRMDecoy::findFixedResidues(sequence, keepN, keepC, keep_const_pattern); }, "sequence"_a, "keepN"_a, "keepC"_a, "keep_const_pattern"_a, 
            R"doc(
Generate decoys from a TargetedExperiment
Will generate decoy peptides for each target peptide provided in exp and
write them into the decoy experiment
Valid methods: shuffle, reverse, pseudo-reverse
If theoretical is true, the target transitions will be returned but their
masses will be adjusted to match the theoretical value of the fragment ion
that is the most likely explanation for the product
`mz_threshold` is used for the matching of theoretical ion series to the observed one
To generate decoys with different precursor mass, use the "switchKR" flag
which switches terminal K/R (switches K to R and R to K). This generates
different precursor m/z and ensures that the y ion series has a different
mass. For a description of the procedure, see (supplemental material)
Bruderer et al. Mol Cell Proteomics. 2017. 10.1074/mcp.RA117.000314.
)doc")
        .def("generateDecoys", [](const OpenMS::MRMDecoy& self,
                                  const OpenMS::TargetedExperiment& exp,
                                  const OpenMS::String& method,
                                  double aim_decoy_fraction,
                                  bool switchKR,
                                  const OpenMS::String& decoy_tag,
                                  int max_attempts,
                                  double identity_threshold,
                                  double precursor_mz_shift,
                                  double product_mz_shift,
                                  double product_mz_threshold,
                                  const std::vector<OpenMS::String>& fragment_types,
                                  const std::vector<size_t>& fragment_charges,
                                  bool enable_specific_losses,
                                  bool enable_unspecific_losses,
                                  int round_decPow) {
            OpenMS::TargetedExperiment dec;
            self.generateDecoys(exp, dec, method, aim_decoy_fraction, switchKR, decoy_tag,
                                max_attempts, identity_threshold, precursor_mz_shift,
                                product_mz_shift, product_mz_threshold, fragment_types,
                                fragment_charges, enable_specific_losses, enable_unspecific_losses,
                                round_decPow);
            return dec;
        }, "exp"_a, "method"_a, "aim_decoy_fraction"_a, "switchKR"_a, "decoy_tag"_a,
           "max_attempts"_a, "identity_threshold"_a, "precursor_mz_shift"_a,
           "product_mz_shift"_a, "product_mz_threshold"_a, "fragment_types"_a,
           "fragment_charges"_a, "enable_specific_losses"_a, "enable_unspecific_losses"_a,
           "round_decPow"_a = -4, "Generate decoys from a TargetedExperiment")
        ;
    def_ProgressLogger<OpenMS::MRMDecoy>(mrmdecoy_class);

    // -----------------------------------------------------------------------
    // MS2File
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MS2File, OpenMS::ProgressLogger>(m, "MS2File", 
        R"doc(
MS2 input file adapter
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("load", [](OpenMS::MS2File& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) { self.load(filename, exp); }, "filename"_a, "exp"_a, "Load MS2 file")
        ;

    // -----------------------------------------------------------------------
    // MSDataCachedConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MSDataCachedConsumer, OpenMS::Internal::CachedMzMLHandler>(m, "MSDataCachedConsumer", 
        R"doc(
Transforming and cached writing consumer of MS data
Is able to transform a spectrum on the fly while it is read using a
function pointer that can be set on the object. The spectra is then
cached to disk using the functions provided in CachedMzMLHandler.
)doc")
        .def(nb::init<OpenMS::String, bool>())
        .def("consumeSpectrum", [](OpenMS::MSDataCachedConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a)
        .def("consumeChromatogram", [](OpenMS::MSDataCachedConsumer& self, OpenMS::MSChromatogram& c) { return self.consumeChromatogram(c); }, "c"_a, 
            R"doc(
Write a spectrum to the output file
May delete data from spectrum (if clearData is set)
)doc")
        .def("setExpectedSize", [](OpenMS::MSDataCachedConsumer& self, size_t p0, size_t p1) { return self.setExpectedSize(p0, p1); })
        .def("setExperimentalSettings", [](OpenMS::MSDataCachedConsumer& self, const OpenMS::ExperimentalSettings& p0) { return self.setExperimentalSettings(p0); }, 
            R"doc(
Write a chromatogram to the output file
May delete data from chromatogram (if clearData is set)
)doc")
        ;

    // -----------------------------------------------------------------------
    // MapAlignmentAlgorithmIdentification
    // -----------------------------------------------------------------------
    auto mapalignmentalgorithmidentification_class = nb::class_<OpenMS::MapAlignmentAlgorithmIdentification, OpenMS::DefaultParamHandler>(m, "MapAlignmentAlgorithmIdentification", 
        R"doc(
A map alignment algorithm based on peptide identifications from MS2
spectra
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("setReference", [](OpenMS::MapAlignmentAlgorithmIdentification& self, const OpenMS::FeatureMap& ref) { self.setReference(ref); }, "ref"_a, "Sets the reference for alignment (FeatureMap)")
        .def("setReference", [](OpenMS::MapAlignmentAlgorithmIdentification& self, const OpenMS::ConsensusMap& ref) { self.setReference(ref); }, "ref"_a, "Sets the reference for alignment (ConsensusMap)")
        .def("align", [](OpenMS::MapAlignmentAlgorithmIdentification& self, const OpenMS::FeatureMap& map) {
            std::vector<OpenMS::FeatureMap> maps = {map};
            std::vector<OpenMS::TransformationDescription> trafos;
            self.align(maps, trafos);
            return trafos.empty() ? OpenMS::TransformationDescription() : trafos[0];
        }, "map"_a, "Aligns a FeatureMap and returns the transformation")
        .def("align", [](OpenMS::MapAlignmentAlgorithmIdentification& self, const OpenMS::ConsensusMap& map) {
            std::vector<OpenMS::ConsensusMap> maps = {map};
            std::vector<OpenMS::TransformationDescription> trafos;
            self.align(maps, trafos);
            return trafos.empty() ? OpenMS::TransformationDescription() : trafos[0];
        }, "map"_a, "Aligns a ConsensusMap and returns the transformation")
        ;
    def_ProgressLogger<OpenMS::MapAlignmentAlgorithmIdentification>(mapalignmentalgorithmidentification_class);

    // -----------------------------------------------------------------------
    // MapAlignmentAlgorithmPoseClustering
    // -----------------------------------------------------------------------
    auto mapalignmentalgorithmposeclustering_class = nb::class_<OpenMS::MapAlignmentAlgorithmPoseClustering, OpenMS::DefaultParamHandler>(m, "MapAlignmentAlgorithmPoseClustering", 
        R"doc(
A map alignment algorithm based on pose clustering
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("setReference", [](OpenMS::MapAlignmentAlgorithmPoseClustering& self, const OpenMS::FeatureMap& ref) { self.setReference(ref); }, "ref"_a, "Sets the reference map for alignment (FeatureMap)")
        .def("setReference", [](OpenMS::MapAlignmentAlgorithmPoseClustering& self, const OpenMS::PeakMap& ref) { self.setReference(ref); }, "ref"_a, "Sets the reference map for alignment (PeakMap)")
        .def("align", [](OpenMS::MapAlignmentAlgorithmPoseClustering& self, const OpenMS::FeatureMap& map) {
            OpenMS::TransformationDescription trafo;
            self.align(map, trafo);
            return trafo;
        }, "map"_a, "Aligns a FeatureMap to the reference and returns the transformation")
        .def("align", [](OpenMS::MapAlignmentAlgorithmPoseClustering& self, const OpenMS::PeakMap& map) {
            OpenMS::TransformationDescription trafo;
            self.align(map, trafo);
            return trafo;
        }, "map"_a, "Aligns a PeakMap to the reference and returns the transformation")
        ;
    def_ProgressLogger<OpenMS::MapAlignmentAlgorithmPoseClustering>(mapalignmentalgorithmposeclustering_class);

    // -----------------------------------------------------------------------
    // MascotGenericFile
    // -----------------------------------------------------------------------
    auto mascotgenericfile_class = nb::class_<OpenMS::MascotGenericFile, OpenMS::ProgressLogger>(m, "MascotGenericFile", 
        R"doc(
Read/write Mascot generic files (MGF)
ProgressLogger
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("store", [](OpenMS::MascotGenericFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& experiment, bool compact) { return self.store(filename, experiment, compact); }, "filename"_a, "experiment"_a, "compact"_a = false)
        .def("store", [](OpenMS::MascotGenericFile& self, std::basic_ostream<char>& os, const OpenMS::String& filename, const OpenMS::MSExperiment& experiment, bool compact) { return self.store(os, filename, experiment, compact); }, "os"_a, "filename"_a, "experiment"_a, "compact"_a = false)
        .def("getHTTPPeakListEnclosure", [](const OpenMS::MascotGenericFile& self, const OpenMS::String& filename) { return self.getHTTPPeakListEnclosure(filename); }, "filename"_a, 
            R"doc(
Loads a Mascot Generic File into a PeakMap
:param filename: File name which the map should be read from
:param exp: The map which is filled with the data from the given file
:raises:
Exception: FileNotFound is thrown if the given file could not be found
)doc")
        .def("load", [](OpenMS::MascotGenericFile& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) { self.load(filename, exp); }, "filename"_a, "exp"_a, "Load MGF file")
        ;
    def_DefaultParamHandler<OpenMS::MascotGenericFile>(mascotgenericfile_class);

    // -----------------------------------------------------------------------
    // MassTraceDetection
    // -----------------------------------------------------------------------
    auto masstracedetection_class = nb::class_<OpenMS::MassTraceDetection, OpenMS::DefaultParamHandler>(m, "MassTraceDetection", 
        R"doc(
A mass trace extraction method that gathers peaks similar in m/z and
moving along retention time
ProgressLogger
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::MassTraceDetection& self, const OpenMS::PeakMap& input, size_t max_traces) {
            std::vector<OpenMS::MassTrace> found_masstraces;
            self.run(input, found_masstraces, max_traces);
            return found_masstraces;
        }, "input"_a, "max_traces"_a = 0,
            R"doc(Main method of MassTraceDetection. Extracts mass traces of a MSExperiment.

:param input: MSExperiment (PeakMap) holding the raw data
:param max_traces: Maximum number of mass traces to extract (0 = no limit)
:return: List of found MassTrace objects
)doc")
        ;
    def_ProgressLogger<OpenMS::MassTraceDetection>(masstracedetection_class);

    // -----------------------------------------------------------------------
    // MasstraceCorrelator
    // -----------------------------------------------------------------------
    auto masstracecorrelator_class = nb::class_<OpenMS::MasstraceCorrelator, OpenMS::DefaultParamHandler>(m, "MasstraceCorrelator", 
        R"doc(
Correlates individual masstraces found in mass spectrometric maps
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("createPseudoSpectra", [](OpenMS::MasstraceCorrelator& self, const OpenMS::ConsensusMap& map, OpenMS::MSExperiment& pseudo_spectra, size_t min_peak_nr, double min_correlation, int max_lag, double max_rt_apex_difference) { return self.createPseudoSpectra(map, pseudo_spectra, min_peak_nr, min_correlation, max_lag, max_rt_apex_difference); }, "map"_a, "pseudo_spectra"_a, "min_peak_nr"_a, "min_correlation"_a, "max_lag"_a, "max_rt_apex_difference"_a)
        ;
    def_ProgressLogger<OpenMS::MasstraceCorrelator>(masstracecorrelator_class);

    // -----------------------------------------------------------------------
    // MetaboliteSpectralMatching
    // -----------------------------------------------------------------------
    auto metabolitespectralmatching_class = nb::class_<OpenMS::MetaboliteSpectralMatching, OpenMS::DefaultParamHandler>(m, "MetaboliteSpectralMatching", "OpenMS class MetaboliteSpectralMatching")
        .def(nb::init<>())
        .def_static("computeHyperScore", [](double fragment_mass_error, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::MSSpectrum& db_spectrum, double mz_lower_bound) { return OpenMS::MetaboliteSpectralMatching::computeHyperScore(fragment_mass_error, fragment_mass_tolerance_unit_ppm, exp_spectrum, db_spectrum, mz_lower_bound); }, "fragment_mass_error"_a, "fragment_mass_tolerance_unit_ppm"_a, "exp_spectrum"_a, "db_spectrum"_a, "mz_lower_bound"_a)
        .def_static("computeHyperScoreWithAnnotations", [](double fragment_mass_error, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::MSSpectrum& db_spectrum, double mz_lower_bound) {
            std::vector<OpenMS::PeptideHit::PeakAnnotation> annotations;
            auto score = OpenMS::MetaboliteSpectralMatching::computeHyperScore(fragment_mass_error, fragment_mass_tolerance_unit_ppm, exp_spectrum, db_spectrum, annotations, mz_lower_bound);
            return nb::make_tuple(score, annotations);
        }, "fragment_mass_error"_a, "fragment_mass_tolerance_unit_ppm"_a, "exp_spectrum"_a, "db_spectrum"_a, "mz_lower_bound"_a)
        .def("run", [](OpenMS::MetaboliteSpectralMatching& self, OpenMS::MSExperiment& p0, OpenMS::MSExperiment& p1, OpenMS::MzTab& p2) {
            OpenMS::String out_spectra;
            self.run(p0, p1, p2, out_spectra);
            return out_spectra;
        })
        ;
    def_ProgressLogger<OpenMS::MetaboliteSpectralMatching>(metabolitespectralmatching_class);

    // -----------------------------------------------------------------------
    // MorphologicalFilter
    // -----------------------------------------------------------------------
    auto morphologicalfilter_class = nb::class_<OpenMS::MorphologicalFilter, OpenMS::ProgressLogger>(m, "MorphologicalFilter", 
        R"doc(
An iterator wrapper to access peak intensities instead of the peak
itself
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("filter", [](OpenMS::MorphologicalFilter& self, OpenMS::MSSpectrum& spectrum) { return self.filter(spectrum); }, "spectrum"_a)
        .def("filterExperiment", [](OpenMS::MorphologicalFilter& self, OpenMS::MSExperiment& exp) { return self.filterExperiment(exp); }, "exp"_a, 
            R"doc(
Applies the morphological filtering operation to an MSSpectrum
If the size of the structuring element is given in 'Thomson', the number of data points for
the structuring element is computed as follows:
- The data points are assumed to be uniformly spaced.  We compute the
average spacing from the position of the first and the last peak and the
total number of peaks in the input range
- The number of data points in the structuring element is computed
from struc_size and the average spacing, and rounded up to an odd
number
)doc")
        ;
    def_DefaultParamHandler<OpenMS::MorphologicalFilter>(morphologicalfilter_class);

    // -----------------------------------------------------------------------
    // OpenPepXLAlgorithm
    // -----------------------------------------------------------------------
    auto openpepxlalgorithm_class = nb::class_<OpenMS::OpenPepXLAlgorithm, OpenMS::DefaultParamHandler>(m, "OpenPepXLAlgorithm", 
        R"doc(
Search for peptide pairs linked with a labeled cross-linker
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::OpenPepXLAlgorithm& self,
                OpenMS::PeakMap& unprocessed_spectra,
                OpenMS::ConsensusMap& cfeatures,
                std::vector<OpenMS::FASTAFile::FASTAEntry>& fasta_db,
                std::vector<OpenMS::ProteinIdentification>& protein_ids,
                OpenMS::PeptideIdentificationList& peptide_ids,
                OpenMS::OPXLDataStructs::PreprocessedPairSpectra& preprocessed_pair_spectra,
                std::vector<std::pair<size_t, size_t>>& spectrum_pairs,
                std::vector<std::vector<OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch>>& all_top_csms,
                OpenMS::PeakMap& spectra) {
            return self.run(unprocessed_spectra, cfeatures, fasta_db, protein_ids, peptide_ids,
                preprocessed_pair_spectra, spectrum_pairs, all_top_csms, spectra);
        }, "unprocessed_spectra"_a, "cfeatures"_a, "fasta_db"_a, "protein_ids"_a,
           "peptide_ids"_a, "preprocessed_pair_spectra"_a, "spectrum_pairs"_a,
           "all_top_csms"_a, "spectra"_a,
           "Run the OpenPepXL cross-link search algorithm")
        ;
    // ExitCodes enum nested under OpenPepXLAlgorithm
    nb::enum_<OpenMS::OpenPepXLAlgorithm::ExitCodes>(openpepxlalgorithm_class, "ExitCodes", nb::is_arithmetic())
        .value("EXECUTION_OK", OpenMS::OpenPepXLAlgorithm::ExitCodes::EXECUTION_OK)
        .value("ILLEGAL_PARAMETERS", OpenMS::OpenPepXLAlgorithm::ExitCodes::ILLEGAL_PARAMETERS)
        .value("UNEXPECTED_RESULT", OpenMS::OpenPepXLAlgorithm::ExitCodes::UNEXPECTED_RESULT)
        .value("INCOMPATIBLE_INPUT_DATA", OpenMS::OpenPepXLAlgorithm::ExitCodes::INCOMPATIBLE_INPUT_DATA)
        ;

    // -----------------------------------------------------------------------
    // PEFFFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PEFFFile, OpenMS::ProgressLogger>(m, "PEFFFile", 
        R"doc(
File adapter for PEFF (PSI Extended FASTA Format) files
PEFF extends FASTA with rich annotations for modifications, variants,
processed regions, and proteoforms.
Usage:
.. code-block:: python
# Batch loading
peff = PEFFFile()
entries = []
headers = []
peff.load("proteins.peff", entries, headers)
for entry in entries:
print(entry.identifier, len(entry.modifications))
# Streaming (memory-efficient for large files)
peff = PEFFFile()
peff.readStart("proteins.peff")
entry = PEFFEntry()
while peff.readNext(entry):
print(entry.identifier)
)doc")
        .def(nb::init<>())
        .def("load", [](const OpenMS::PEFFFile& self, const OpenMS::String& filename) { std::vector<OpenMS::PEFFEntry> entries; std::vector<OpenMS::PEFFDatabaseMetadata> headers; self.load(filename, entries, headers); return std::make_tuple(entries, headers); }, "filename"_a)
        .def("readStart", [](OpenMS::PEFFFile& self, const OpenMS::String& filename) { return self.readStart(filename); }, "filename"_a, 
            R"doc(
Stores entries to a PEFF file with the given header
:param filename: The output file path
:param entries: The entries to store
:param header: The database metadata header
)doc")
        .def("readNext", [](OpenMS::PEFFFile& self, OpenMS::PEFFEntry& entry) { return self.readNext(entry); }, "entry"_a, 
            R"doc(
Prepares a PEFF file for streamed reading using readNext()
:raises:
Exception:FileNotFound is thrown if the file does not exist
Exception:FileNotReadable is thrown if the file cannot be read
)doc")
        .def("getHeaders", [](const OpenMS::PEFFFile& self) -> const std::vector<OpenMS::PEFFDatabaseMetadata> & { return self.getHeaders(); }, nb::rv_policy::reference_internal, 
            R"doc(
Reads the next PEFF entry from file
:return: True if entry was read; False if EOF was reached
)doc")
        .def("atEnd", [](const OpenMS::PEFFFile& self) { return self.atEnd(); }, "Returns the headers parsed during readStart()")
        .def("writeStart", [](OpenMS::PEFFFile& self, const OpenMS::String& filename, const OpenMS::PEFFDatabaseMetadata& header) { return self.writeStart(filename, header); }, "filename"_a, "header"_a, "Returns True if the end of the file has been reached")
        .def("writeStart", [](OpenMS::PEFFFile& self, const OpenMS::String& filename, const std::vector<OpenMS::PEFFDatabaseMetadata>& headers) { return self.writeStart(filename, headers); }, "filename"_a, "headers"_a, "Returns True if the end of the file has been reached")
        .def("writeNext", [](OpenMS::PEFFFile& self, const OpenMS::PEFFEntry& entry) { return self.writeNext(entry); }, "entry"_a, 
            R"doc(
Prepares a PEFF file for streamed writing using writeNext()
:raises:
Exception:UnableToCreateFile is thrown if the file cannot be created
)doc")
        .def("writeEnd", [](OpenMS::PEFFFile& self) { return self.writeEnd(); }, "Writes the next PEFF entry to the file")
        .def_static("isPEFFFile", [](const OpenMS::String& filename) { return OpenMS::PEFFFile::isPEFFFile(filename); }, "filename"_a)
        .def_static("toProForma", [](const OpenMS::PEFFEntry& entry) { return OpenMS::PEFFFile::toProForma(entry); }, "entry"_a)
        .def("store", [](const OpenMS::PEFFFile& self, const OpenMS::String& filename, const std::vector<OpenMS::PEFFEntry>& entries, const OpenMS::PEFFDatabaseMetadata& header) { self.store(filename, entries, header); }, "filename"_a, "entries"_a, "header"_a, "Store PEFF file with single header")
        .def("store", [](const OpenMS::PEFFFile& self, const OpenMS::String& filename, const std::vector<OpenMS::PEFFEntry>& entries, const std::vector<OpenMS::PEFFDatabaseMetadata>& headers) { self.store(filename, entries, headers); }, "filename"_a, "entries"_a, "headers"_a, "Store PEFF file with multiple headers")
        ;

    // -----------------------------------------------------------------------
    // PeakBoundary (PeakPickerHiRes::PeakBoundary)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakPickerHiRes::PeakBoundary>(m, "PeakBoundary",
        "Peak boundary information from PeakPickerHiRes")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::PeakPickerHiRes::PeakBoundary& self) { return OpenMS::PeakPickerHiRes::PeakBoundary(self); })
        .def("__deepcopy__", [](const OpenMS::PeakPickerHiRes::PeakBoundary& self, nb::dict) { return OpenMS::PeakPickerHiRes::PeakBoundary(self); }, "memo"_a)
        .def_rw("mz_min", &OpenMS::PeakPickerHiRes::PeakBoundary::mz_min)
        .def_rw("mz_max", &OpenMS::PeakPickerHiRes::PeakBoundary::mz_max)
        ;

    // -----------------------------------------------------------------------
    // PeakPickerHiRes
    // -----------------------------------------------------------------------
    auto peakpickerhires_class = nb::class_<OpenMS::PeakPickerHiRes, OpenMS::DefaultParamHandler>(m, "PeakPickerHiRes", 
        R"doc(
This class implements a fast peak-picking algorithm best suited for
high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution
data, the signals of ions with similar mass-to-charge ratios (m/z)
exhibit little or no overlapping and therefore allow for a clear
separation. Furthermore, ion signals tend to show well-defined peak
shapes with narrow peak width
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("pick", [](const OpenMS::PeakPickerHiRes& self, const OpenMS::MSSpectrum& input, OpenMS::MSSpectrum& output) { return self.pick(input, output); }, "input"_a, "output"_a)
        .def("pick", [](const OpenMS::PeakPickerHiRes& self, const OpenMS::MSChromatogram& input, OpenMS::MSChromatogram& output) { return self.pick(input, output); }, "input"_a, "output"_a)
        .def("pick", [](const OpenMS::PeakPickerHiRes& self, const OpenMS::Mobilogram& input, OpenMS::Mobilogram& output) { return self.pick(input, output); }, "input"_a, "output"_a)
        .def("pick", [](const OpenMS::PeakPickerHiRes& self, const OpenMS::MSSpectrum& input, OpenMS::MSSpectrum& output, bool check_spacings) {
            std::vector<OpenMS::PeakPickerHiRes::PeakBoundary> boundaries;
            self.pick(input, output, boundaries, check_spacings);
            return boundaries;
        }, "input"_a, "output"_a, "check_spacings"_a)
        .def("pick", [](const OpenMS::PeakPickerHiRes& self, const OpenMS::MSChromatogram& input, OpenMS::MSChromatogram& output, bool check_spacings) {
            std::vector<OpenMS::PeakPickerHiRes::PeakBoundary> boundaries;
            self.pick(input, output, boundaries, check_spacings);
            return boundaries;
        }, "input"_a, "output"_a, "check_spacings"_a)
        .def("pick", [](const OpenMS::PeakPickerHiRes& self, const OpenMS::Mobilogram& input, OpenMS::Mobilogram& output, bool check_spacings) {
            std::vector<OpenMS::PeakPickerHiRes::PeakBoundary> boundaries;
            self.pick(input, output, boundaries, check_spacings);
            return boundaries;
        }, "input"_a, "output"_a, "check_spacings"_a)
        .def("pickExperiment", [](const OpenMS::PeakPickerHiRes& self, OpenMS::PeakMap& input, bool check_spectrum_type) {
            OpenMS::PeakMap output;
            self.pickExperiment(input, output, check_spectrum_type);
            return output;
        }, "input"_a, "check_spectrum_type"_a = true, "Picks peaks in a full experiment (MSExperiment)")
        ;
    def_ProgressLogger<OpenMS::PeakPickerHiRes>(peakpickerhires_class);

    // -----------------------------------------------------------------------
    // PeakPickerIterative
    // -----------------------------------------------------------------------
    auto peakpickeriterative_class = nb::class_<OpenMS::PeakPickerIterative, OpenMS::DefaultParamHandler>(m, "PeakPickerIterative", 
        R"doc(
Iterative peak picker that uses seed-based centroiding to detect and
integrate peaks in profile spectra
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("pick", [](OpenMS::PeakPickerIterative& self, const OpenMS::MSSpectrum& input, OpenMS::MSSpectrum& output) { return self.pick(input, output); }, "input"_a, "output"_a)
        .def("pickExperiment", [](OpenMS::PeakPickerIterative& self, const OpenMS::MSExperiment& input, OpenMS::MSExperiment& output) { return self.pickExperiment(input, output); }, "input"_a, "output"_a, 
            R"doc(
This will pick one single spectrum. The PeakPickerHiRes is used to
generate seeds, these seeds are then used to re-center the mass and
compute peak width and integrated intensity of the peak
Finally, other peaks that would fall within the primary peak are
discarded
The output are the remaining peaks
)doc")
        ;
    def_ProgressLogger<OpenMS::PeakPickerIterative>(peakpickeriterative_class);

    // -----------------------------------------------------------------------
    // PeptideIndexing
    // -----------------------------------------------------------------------
    auto peptideindexing_class = nb::class_<OpenMS::PeptideIndexing, OpenMS::DefaultParamHandler>(m, "PeptideIndexing", 
        R"doc(
DefaultParamHandler

Refreshes the protein references for all peptide hits in a vector of PeptideIdentifications and adds target/decoy information
All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". For proteins the possible values are "target" and "decoy",
depending on whether the protein accession contains the decoy pattern (parameter `decoy_string`) as a suffix or prefix, respectively (see parameter `prefix`).
For peptides, the possible values are "target", "decoy" and "target+decoy", depending on whether the peptide sequence is found only in target proteins,
only in decoy proteins, or in both. The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool.
(For FDR calculations, "target+decoy" peptide hits count as target hits.)
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::PeptideIndexing& self, std::vector<OpenMS::FASTAFile::FASTAEntry> proteins, std::vector<OpenMS::ProteinIdentification> prot_ids, OpenMS::PeptideIdentificationList& pep_ids) {
            decltype(self.run(proteins, prot_ids, pep_ids)) result;
            { nb::gil_scoped_release release; result = self.run(proteins, prot_ids, pep_ids); }
            return nb::make_tuple(result, proteins, prot_ids);
        }, "proteins"_a, "prot_ids"_a, "pep_ids"_a)
        .def("run", [](OpenMS::PeptideIndexing& self, OpenMS::FASTAContainer<OpenMS::TFI_File>& proteins, std::vector<OpenMS::ProteinIdentification> prot_ids, OpenMS::PeptideIdentificationList& pep_ids) {
            decltype(self.run(proteins, prot_ids, pep_ids)) result;
            { nb::gil_scoped_release release; result = self.run(proteins, prot_ids, pep_ids); }
            return nb::make_tuple(result, prot_ids);
        }, "proteins"_a, "prot_ids"_a, "pep_ids"_a)
        .def("run", [](OpenMS::PeptideIndexing& self, OpenMS::FASTAContainer<OpenMS::TFI_Vector>& proteins, std::vector<OpenMS::ProteinIdentification> prot_ids, OpenMS::PeptideIdentificationList& pep_ids) {
            decltype(self.run(proteins, prot_ids, pep_ids)) result;
            { nb::gil_scoped_release release; result = self.run(proteins, prot_ids, pep_ids); }
            return nb::make_tuple(result, prot_ids);
        }, "proteins"_a, "prot_ids"_a, "pep_ids"_a)
        .def("getDecoyString", [](const OpenMS::PeptideIndexing& self) { return self.getDecoyString(); })
        .def("isPrefix", [](const OpenMS::PeptideIndexing& self) { return self.isPrefix(); })
        ;
    // ExitCodes enum nested under PeptideIndexing
    nb::enum_<OpenMS::PeptideIndexing::ExitCodes>(peptideindexing_class, "ExitCodes", nb::is_arithmetic())
        .value("EXECUTION_OK", OpenMS::PeptideIndexing::ExitCodes::EXECUTION_OK)
        .value("DATABASE_EMPTY", OpenMS::PeptideIndexing::ExitCodes::DATABASE_EMPTY)
        .value("PEPTIDE_IDS_EMPTY", OpenMS::PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY)
        .value("ILLEGAL_PARAMETERS", OpenMS::PeptideIndexing::ExitCodes::ILLEGAL_PARAMETERS)
        .value("UNEXPECTED_RESULT", OpenMS::PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
        ;

    // -----------------------------------------------------------------------
    // OpenSearchModificationAnalysis — nested structs as top-level classes
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSearchModificationAnalysis::ModificationPattern>(m, "ModificationPattern",
        "Stores details of a modification pattern found in the data")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis::ModificationPattern&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis::ModificationPattern& self) { return OpenMS::OpenSearchModificationAnalysis::ModificationPattern(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis::ModificationPattern& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis::ModificationPattern(self); }, "memo"_a)
        .def_rw("count", &OpenMS::OpenSearchModificationAnalysis::ModificationPattern::count)
        .def_rw("masses", &OpenMS::OpenSearchModificationAnalysis::ModificationPattern::masses)
        .def_rw("num_charge_states", &OpenMS::OpenSearchModificationAnalysis::ModificationPattern::num_charge_states)
        ;

    nb::class_<OpenMS::OpenSearchModificationAnalysis::ModificationSummary>(m, "ModificationSummary",
        "Modification summary output with count, name, charge states, and masses")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis::ModificationSummary&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis::ModificationSummary& self) { return OpenMS::OpenSearchModificationAnalysis::ModificationSummary(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis::ModificationSummary& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis::ModificationSummary(self); }, "memo"_a)
        .def_rw("count", &OpenMS::OpenSearchModificationAnalysis::ModificationSummary::count)
        .def_rw("name", &OpenMS::OpenSearchModificationAnalysis::ModificationSummary::name)
        .def_rw("num_charge_states", &OpenMS::OpenSearchModificationAnalysis::ModificationSummary::num_charge_states)
        .def_rw("masses", &OpenMS::OpenSearchModificationAnalysis::ModificationSummary::masses)
        ;

    nb::class_<OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry>(m, "DeltaMassEntry",
        "Statistics for a single delta mass bin in the histogram")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry& self) { return OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry(self); }, "memo"_a)
        .def_rw("delta_mass", &OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry::delta_mass)
        .def_rw("count", &OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry::count)
        .def_rw("unique_peptides", &OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry::unique_peptides)
        .def_rw("num_charge_states", &OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry::num_charge_states)
        .def_rw("percentage", &OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry::percentage)
        .def_rw("mapped_modification", &OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry::mapped_modification)
        .def_rw("is_known_modification", &OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry::is_known_modification)
        ;

    nb::class_<OpenMS::OpenSearchModificationAnalysis::PTMEntry>(m, "PTMEntry",
        "Statistics for a mapped PTM with residue localization")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis::PTMEntry&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis::PTMEntry& self) { return OpenMS::OpenSearchModificationAnalysis::PTMEntry(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis::PTMEntry& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis::PTMEntry(self); }, "memo"_a)
        .def_rw("name", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::name)
        .def_rw("theoretical_mass", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::theoretical_mass)
        .def_rw("observed_mass", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::observed_mass)
        .def_rw("mass_deviation", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::mass_deviation)
        .def_rw("count", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::count)
        .def_rw("unique_peptides", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::unique_peptides)
        .def_rw("num_charge_states", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::num_charge_states)
        .def_rw("percentage", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::percentage)
        .def_rw("residue_counts", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::residue_counts)
        .def_rw("target_residues", &OpenMS::OpenSearchModificationAnalysis::PTMEntry::target_residues)
        ;

    nb::class_<OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics>(m, "DeltaMassStatistics",
        "Container for delta mass statistics table")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics& self) { return OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics(self); }, "memo"_a)
        .def_rw("entries", &OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics::entries)
        .def_rw("total_psms", &OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics::total_psms)
        .def_rw("modified_psms", &OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics::modified_psms)
        .def_rw("unmodified_psms", &OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics::unmodified_psms)
        .def_rw("mean_delta_mass", &OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics::mean_delta_mass)
        .def_rw("median_delta_mass", &OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics::median_delta_mass)
        ;

    nb::class_<OpenMS::OpenSearchModificationAnalysis::PTMStatistics>(m, "PTMStatistics",
        "Container for PTM statistics table")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis::PTMStatistics&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis::PTMStatistics& self) { return OpenMS::OpenSearchModificationAnalysis::PTMStatistics(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis::PTMStatistics& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis::PTMStatistics(self); }, "memo"_a)
        .def_rw("entries", &OpenMS::OpenSearchModificationAnalysis::PTMStatistics::entries)
        .def_rw("total_modified_psms", &OpenMS::OpenSearchModificationAnalysis::PTMStatistics::total_modified_psms)
        .def_rw("unknown_modification_psms", &OpenMS::OpenSearchModificationAnalysis::PTMStatistics::unknown_modification_psms)
        .def_rw("num_unique_modifications", &OpenMS::OpenSearchModificationAnalysis::PTMStatistics::num_unique_modifications)
        ;

    nb::class_<OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult>(m, "OpenSearchAnalysisResult",
        "Combined result of open search modification analysis")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult& self) { return OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult(self); }, "memo"_a)
        .def_rw("delta_mass_stats", &OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult::delta_mass_stats)
        .def_rw("ptm_stats", &OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult::ptm_stats)
        .def_rw("summaries", &OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult::summaries)
        ;

    // -----------------------------------------------------------------------
    // OpenSearchModificationAnalysis
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSearchModificationAnalysis>(m, "OpenSearchModificationAnalysis",
        "Utility class for analyzing modification patterns in open search results")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSearchModificationAnalysis&>())
        .def("__copy__", [](const OpenMS::OpenSearchModificationAnalysis& self) { return OpenMS::OpenSearchModificationAnalysis(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSearchModificationAnalysis& self, nb::dict) { return OpenMS::OpenSearchModificationAnalysis(self); }, "memo"_a)
        .def("analyzeModifications", &OpenMS::OpenSearchModificationAnalysis::analyzeModifications,
            "peptide_ids"_a, "precursor_mass_tolerance"_a = 5.0,
            "precursor_mass_tolerance_unit_ppm"_a = true, "use_smoothing"_a = false,
            "output_file"_a = OpenMS::String(""),
            "Complete analysis workflow: analyze patterns and map to modifications")
        .def("analyzeModificationsWithStatistics", &OpenMS::OpenSearchModificationAnalysis::analyzeModificationsWithStatistics,
            "peptide_ids"_a, "precursor_mass_tolerance"_a = 5.0,
            "precursor_mass_tolerance_unit_ppm"_a = true, "use_smoothing"_a = false,
            "output_file"_a = OpenMS::String(""),
            "Complete analysis returning structured statistics tables")
        .def("generatePTMStatistics", &OpenMS::OpenSearchModificationAnalysis::generatePTMStatistics,
            "peptide_ids"_a, "precursor_mass_tolerance"_a = 5.0,
            "precursor_mass_tolerance_unit_ppm"_a = true,
            "Generate PTM statistics table with residue localization")
        .def("analyzeResidueFrequency", &OpenMS::OpenSearchModificationAnalysis::analyzeResidueFrequency,
            "peptide_ids"_a, "delta_mass"_a, "tolerance"_a = 0.01,
            "Analyze which amino acid residues are associated with a delta mass")
        .def("writeDeltaMassStatistics", &OpenMS::OpenSearchModificationAnalysis::writeDeltaMassStatistics,
            "stats"_a, "output_file"_a,
            "Write delta mass statistics to a TSV file")
        .def("writePTMStatistics", &OpenMS::OpenSearchModificationAnalysis::writePTMStatistics,
            "stats"_a, "output_file"_a,
            "Write PTM statistics to a TSV file")
        ;

    // -----------------------------------------------------------------------
    // ProSEAlgorithm
    // -----------------------------------------------------------------------
    // SearchResult struct (nested in ProSEAlgorithm)
    nb::class_<OpenMS::ProSEAlgorithm::SearchResult>(m, "SearchResult",
        "Comprehensive search result including modification analysis")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProSEAlgorithm::SearchResult&>())
        .def("__copy__", [](const OpenMS::ProSEAlgorithm::SearchResult& self) { return OpenMS::ProSEAlgorithm::SearchResult(self); })
        .def("__deepcopy__", [](const OpenMS::ProSEAlgorithm::SearchResult& self, nb::dict) { return OpenMS::ProSEAlgorithm::SearchResult(self); }, "memo"_a)
        .def_rw("exit_code", &OpenMS::ProSEAlgorithm::SearchResult::exit_code)
        .def_rw("protein_ids", &OpenMS::ProSEAlgorithm::SearchResult::protein_ids)
        .def_rw("peptide_ids", &OpenMS::ProSEAlgorithm::SearchResult::peptide_ids)
        .def_rw("modification_analysis", &OpenMS::ProSEAlgorithm::SearchResult::modification_analysis)
        .def_rw("is_open_search", &OpenMS::ProSEAlgorithm::SearchResult::is_open_search)
        ;

    auto prosealgorithm_class = nb::class_<OpenMS::ProSEAlgorithm, OpenMS::DefaultParamHandler>(m, "ProSEAlgorithm",
        R"doc(
DefaultParamHandler
ProgressLogger

Fragment-index-based peptide database search algorithm (experimental).
Provides a self-contained search engine that matches MS/MS spectra against a protein
database using an FI (Fragment Index). Typical usage:
- Configure parameters via DefaultParamHandler (mass tolerances, enzyme, charges, etc.)
- Call search() with an input mzML file and a FASTA database to populate identification
outputs (ProteinIdentification and PeptideIdentificationList)
- Intended for educational/prototyping use and to demonstrate FI-backed searching
)doc")
        .def(nb::init<>())
        // in-memory search with prot_ids output parameter (4-arg, most specific — must be first)
        .def("search", [](const OpenMS::ProSEAlgorithm& self, OpenMS::PeakMap& spectra, const std::vector<OpenMS::FASTAFile::FASTAEntry>& fasta_db, nb::list prot_ids_out, OpenMS::PeptideIdentificationList& pep_ids) {
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            auto result = self.search(spectra, fasta_db, prot_ids, pep_ids);
            for (auto& p : prot_ids) {
                prot_ids_out.append(nb::cast(std::move(p)));
            }
            return result;
        }, "spectra"_a, "fasta_db"_a, "prot_ids"_a, "pep_ids"_a,
           "In-memory search. Populates prot_ids list and returns ExitCodes")
        // in-memory search overload (PeakMap + FASTAEntry vector, 3-arg returns tuple)
        .def("search", [](const OpenMS::ProSEAlgorithm& self, OpenMS::PeakMap& spectra, const std::vector<OpenMS::FASTAFile::FASTAEntry>& fasta_db, OpenMS::PeptideIdentificationList& pep_ids) {
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            auto result = self.search(spectra, fasta_db, prot_ids, pep_ids);
            return nb::make_tuple(result, prot_ids);
        }, "spectra"_a, "fasta_db"_a, "pep_ids"_a,
           "In-memory search. Returns (ExitCodes, list[ProteinIdentification])")
        // file-based search with prot_ids output parameter (4-arg)
        .def("search", [](const OpenMS::ProSEAlgorithm& self, const OpenMS::String& in_mzML, const OpenMS::String& in_db, nb::list prot_ids_out, OpenMS::PeptideIdentificationList& pep_ids) {
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            auto result = self.search(in_mzML, in_db, prot_ids, pep_ids);
            for (auto& p : prot_ids) {
                prot_ids_out.append(nb::cast(std::move(p)));
            }
            return result;
        }, "in_mzML"_a, "in_db"_a, "prot_ids"_a, "pep_ids"_a,
           "File-based search. Populates prot_ids list and returns ExitCodes")
        // file-based search (3-arg returns tuple)
        .def("search", [](const OpenMS::ProSEAlgorithm& self, const OpenMS::String& in_mzML, const OpenMS::String& in_db, OpenMS::PeptideIdentificationList& pep_ids) {
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            auto result = self.search(in_mzML, in_db, prot_ids, pep_ids);
            return nb::make_tuple(result, prot_ids);
        }, "in_mzML"_a, "in_db"_a, "pep_ids"_a,
           "File-based search. Returns (ExitCodes, list[ProteinIdentification])")
        // in-memory searchWithModificationAnalysis (more specific types — must be first)
        .def("searchWithModificationAnalysis",
            static_cast<OpenMS::ProSEAlgorithm::SearchResult (OpenMS::ProSEAlgorithm::*)(OpenMS::PeakMap&, const std::vector<OpenMS::FASTAFile::FASTAEntry>&, const OpenMS::String&) const>(
                &OpenMS::ProSEAlgorithm::searchWithModificationAnalysis),
            "spectra"_a, "fasta_db"_a, "output_base_name"_a = OpenMS::String(""),
            "In-memory search with modification analysis. Returns SearchResult")
        // file-based searchWithModificationAnalysis
        .def("searchWithModificationAnalysis",
            static_cast<OpenMS::ProSEAlgorithm::SearchResult (OpenMS::ProSEAlgorithm::*)(const OpenMS::String&, const OpenMS::String&, const OpenMS::String&) const>(
                &OpenMS::ProSEAlgorithm::searchWithModificationAnalysis),
            "in_mzML"_a, "in_db"_a, "output_base_name"_a = OpenMS::String(""),
            "File-based search with modification analysis. Returns SearchResult")
        ;
    def_ProgressLogger<OpenMS::ProSEAlgorithm>(prosealgorithm_class);
    // ProSEAlgorithm_ExitCodes enum nested under ProSEAlgorithm
    nb::enum_<OpenMS::ProSEAlgorithm::ExitCodes>(prosealgorithm_class, "ProSEAlgorithm_ExitCodes", nb::is_arithmetic())
        .value("EXECUTION_OK", OpenMS::ProSEAlgorithm::ExitCodes::EXECUTION_OK)
        .value("INPUT_FILE_EMPTY", OpenMS::ProSEAlgorithm::ExitCodes::INPUT_FILE_EMPTY)
        .value("UNEXPECTED_RESULT", OpenMS::ProSEAlgorithm::ExitCodes::UNEXPECTED_RESULT)
        .value("UNKNOWN_ERROR", OpenMS::ProSEAlgorithm::ExitCodes::UNKNOWN_ERROR)
        .value("ILLEGAL_PARAMETERS", OpenMS::ProSEAlgorithm::ExitCodes::ILLEGAL_PARAMETERS)
        .export_values();

    // -----------------------------------------------------------------------
    // QTClusterFinder
    // -----------------------------------------------------------------------
    auto qtclusterfinder_class = nb::class_<OpenMS::QTClusterFinder, OpenMS::BaseGroupFinder>(m, "QTClusterFinder", 
        R"doc(
A variant of QT clustering for the detection of feature groups
BaseGroupFinder
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::QTClusterFinder& self, const std::vector<OpenMS::FeatureMap>& input_maps, OpenMS::ConsensusMap& result_map) { return self.run(input_maps, result_map); }, "input_maps"_a, "result_map"_a)
        ;
    def_ProgressLogger<OpenMS::QTClusterFinder>(qtclusterfinder_class);

    // -----------------------------------------------------------------------
    // RankScaler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RankScaler, OpenMS::DefaultParamHandler>(m, "RankScaler", 
        R"doc(
Scales each peak by ranking the peaks per spectrum and assigning
intensity according to rank
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RankScaler &>())
        .def("__copy__", [](const OpenMS::RankScaler& self) { return OpenMS::RankScaler(self); })
        .def("__deepcopy__", [](const OpenMS::RankScaler& self, nb::dict) { return OpenMS::RankScaler(self); }, "memo"_a)
        .def("filterPeakSpectrum", [](OpenMS::RankScaler& self, OpenMS::MSSpectrum& spectrum) { return self.filterPeakSpectrum(spectrum); }, "spectrum"_a)
        .def("filterPeakMap", [](OpenMS::RankScaler& self, OpenMS::MSExperiment& exp) { return self.filterPeakMap(exp); }, "exp"_a)
        ;

    // -----------------------------------------------------------------------
    // SavitzkyGolayFilter
    // -----------------------------------------------------------------------
    auto savitzkygolayfilter_class = nb::class_<OpenMS::SavitzkyGolayFilter, OpenMS::ProgressLogger>(m, "SavitzkyGolayFilter", 
        R"doc(
Computes the Savitzky-Golay filter coefficients using QR
decomposition
DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("filter", [](OpenMS::SavitzkyGolayFilter& self, OpenMS::MSSpectrum& spectrum) { return self.filter(spectrum); }, "spectrum"_a, "Removed the noise from an MSSpectrum containing profile data")
        .def("filter", [](OpenMS::SavitzkyGolayFilter& self, OpenMS::MSChromatogram& chromatogram) { return self.filter(chromatogram); }, "chromatogram"_a, "Removed the noise from an MSSpectrum containing profile data")
        .def("filter", [](OpenMS::SavitzkyGolayFilter& self, OpenMS::Mobilogram& mobilogram) { return self.filter(mobilogram); }, "mobilogram"_a, "Removed the noise from an MSSpectrum containing profile data")
        .def("filterExperiment", [](OpenMS::SavitzkyGolayFilter& self, OpenMS::MSExperiment& map) { return self.filterExperiment(map); }, "map"_a, "Removed the noise from an MSExperiment containing profile data")
        ;
    def_DefaultParamHandler<OpenMS::SavitzkyGolayFilter>(savitzkygolayfilter_class);

    // -----------------------------------------------------------------------
    // SimpleSearchEngineAlgorithm
    // -----------------------------------------------------------------------
    auto simplesearchenginealgorithm_class = nb::class_<OpenMS::SimpleSearchEngineAlgorithm, OpenMS::DefaultParamHandler>(m, "SimpleSearchEngineAlgorithm", "OpenMS class SimpleSearchEngineAlgorithm")
        .def(nb::init<>())
        .def("search", [](const OpenMS::SimpleSearchEngineAlgorithm& self, const OpenMS::String& in_mzML, const OpenMS::String& in_db, OpenMS::PeptideIdentificationList& pep_ids) {
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            decltype(self.search(in_mzML, in_db, prot_ids, pep_ids)) result;
            { nb::gil_scoped_release release; result = self.search(in_mzML, in_db, prot_ids, pep_ids); }
            return nb::make_tuple(result, prot_ids);
        }, "in_mzML"_a, "in_db"_a, "pep_ids"_a)
        ;
    def_ProgressLogger<OpenMS::SimpleSearchEngineAlgorithm>(simplesearchenginealgorithm_class);

    // -----------------------------------------------------------------------
    // SimpleTSGXLMS
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SimpleTSGXLMS, OpenMS::DefaultParamHandler>(m, "SimpleTSGXLMS", 
        R"doc(
Generates theoretical spectra for cross-linked peptides
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SimpleTSGXLMS &>())
        .def("__copy__", [](const OpenMS::SimpleTSGXLMS& self) { return OpenMS::SimpleTSGXLMS(self); })
        .def("__deepcopy__", [](const OpenMS::SimpleTSGXLMS& self, nb::dict) { return OpenMS::SimpleTSGXLMS(self); }, "memo"_a)
        .def("getLinearIonSpectrum", [](const OpenMS::SimpleTSGXLMS& self, size_t link_pos, int charge, size_t link_pos_2) { std::vector<OpenMS::SimpleTSGXLMS::SimplePeak> spectrum; OpenMS::AASequence peptide; self.getLinearIonSpectrum(spectrum, peptide, link_pos, charge, link_pos_2); return std::make_tuple(spectrum, peptide); }, "link_pos"_a, "charge"_a, "link_pos_2"_a)
        .def("getXLinkIonSpectrum", [](const OpenMS::SimpleTSGXLMS& self, size_t link_pos, double precursor_mass, int min_charge, int max_charge, size_t link_pos_2) { std::vector<OpenMS::SimpleTSGXLMS::SimplePeak> spectrum; OpenMS::AASequence peptide; self.getXLinkIonSpectrum(spectrum, peptide, link_pos, precursor_mass, min_charge, max_charge, link_pos_2); return std::make_tuple(spectrum, peptide); }, "link_pos"_a, "precursor_mass"_a, "min_charge"_a, "max_charge"_a, "link_pos_2"_a, 
            R"doc(
Generates fragment ions not containing the cross-linker for one peptide
B-ions are generated from the beginning of the peptide up to the first linked position,
y-ions are generated from the second linked position up the end of the peptide
If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
The generated ion types and other additional settings are determined by the tool parameters
:param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
:param peptide: The peptide to fragment
:param link_pos: The position of the cross-linker on the given peptide
:param charge: The maximal charge of the ions
:param link_pos_2: A second position for the linker, in case it is a loop link
)doc")
        .def("getXLinkIonSpectrum", [](const OpenMS::SimpleTSGXLMS& self, bool frag_alpha, int min_charge, int max_charge) { std::vector<OpenMS::SimpleTSGXLMS::SimplePeak> spectrum; OpenMS::OPXLDataStructs::ProteinProteinCrossLink crosslink; self.getXLinkIonSpectrum(spectrum, crosslink, frag_alpha, min_charge, max_charge); return std::make_tuple(spectrum, crosslink); }, "frag_alpha"_a, "min_charge"_a, "max_charge"_a, 
            R"doc(
Generates fragment ions containing the cross-linker for one peptide
B-ions are generated from the first linked position up to the end of the peptide,
y-ions are generated from the beginning of the peptide up to the second linked position
If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
a precursor mass for the two peptides and the linker is needed
In the case of a loop link the precursor mass is the mass of the only peptide and the linker
Although this function is more general, currently it is mainly used for loop-links and mono-links,
because residues in the second, unknown peptide cannot be considered for possible neutral losses
The generated ion types and other additional settings are determined by the tool parameters
:param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
:param peptide: The peptide to fragment
:param link_pos: The position of the cross-linker on the given peptide
:param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
:param mincharge: The minimal charge of the ions
:param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
:param link_pos_2: A second position for the linker, in case it is a loop link
)doc")
        ;

    // -----------------------------------------------------------------------
    // SiriusExportAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SiriusExportAlgorithm, OpenMS::DefaultParamHandler>(m, "SiriusExportAlgorithm", "OpenMS class SiriusExportAlgorithm")
        .def(nb::init<>())
        .def("isFeatureOnly", [](const OpenMS::SiriusExportAlgorithm& self) { return self.isFeatureOnly(); })
        .def("getFilterByNumMassTraces", [](const OpenMS::SiriusExportAlgorithm& self) { return self.getFilterByNumMassTraces(); })
        .def("getPrecursorMzTolerance", [](const OpenMS::SiriusExportAlgorithm& self) { return self.getPrecursorMzTolerance(); })
        .def("getPrecursorRtTolerance", [](const OpenMS::SiriusExportAlgorithm& self) { return self.getPrecursorRtTolerance(); })
        .def("precursorMzToleranceUnitIsPPM", [](const OpenMS::SiriusExportAlgorithm& self) { return self.precursorMzToleranceUnitIsPPM(); })
        .def("isNoMasstraceInfoIsotopePattern", [](const OpenMS::SiriusExportAlgorithm& self) { return self.isNoMasstraceInfoIsotopePattern(); })
        .def("getIsotopePatternIterations", [](const OpenMS::SiriusExportAlgorithm& self) { return self.getIsotopePatternIterations(); })
        .def("preprocessing", [](const OpenMS::SiriusExportAlgorithm& self, const OpenMS::String& featureinfo, const OpenMS::MSExperiment& spectra, OpenMS::FeatureMapping::FeatureMappingInfo& feature_mapping_info, OpenMS::FeatureMapping::FeatureToMs2Indices& feature_ms2_indices) { return self.preprocessing(featureinfo, spectra, feature_mapping_info, feature_ms2_indices); }, "featureinfo"_a, "spectra"_a, "feature_mapping_info"_a, "feature_ms2_indices"_a)
        .def("logFeatureSpectraNumber", [](const OpenMS::SiriusExportAlgorithm& self, const OpenMS::String& featureinfo, const OpenMS::FeatureMapping::FeatureToMs2Indices& feature_ms2_indices, const OpenMS::MSExperiment& spectra) { return self.logFeatureSpectraNumber(featureinfo, feature_ms2_indices, spectra); }, "featureinfo"_a, "feature_ms2_indices"_a, "spectra"_a, 
            R"doc(
Preprocessing needed for SIRIUS
Filter number of masstraces and perform feature mapping
:param featureXML_path: Path to featureXML
:param spectra: Input of MSExperiment with spectra information
:param feature_mapping_info: Emtpy - stores FeatureMaps and KDTreeMaps internally
:param feature_ms2_indices: Empty FeatureToMs2Indices
)doc")
        .def("run", [](const OpenMS::SiriusExportAlgorithm& self, const std::vector<OpenMS::String>& mzML_files, const std::vector<OpenMS::String>& featureXML_files, const OpenMS::String& out_ms, const OpenMS::String& out_compoundinfo) { return self.run(mzML_files, featureXML_files, out_ms, out_compoundinfo); }, "mzML_files"_a, "featureXML_files"_a, "out_ms"_a, "out_compoundinfo"_a, 
            R"doc(
Logs number of features and spectra used
Prints the number of features and spectra used (OPENMS_LOG_INFO)
:param featureXML_path: Path to featureXML
:param feature_ms2_indices: FeatureToMs2Indices with feature mapping
:param spectra: Input of MSExperiment with spectra information
)doc")
        ;

    // -----------------------------------------------------------------------
    // SpectraMerger
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectraMerger, OpenMS::DefaultParamHandler>(m, "SpectraMerger", 
        R"doc(
Offers spectra merging and averaging algorithms to increase the
quality of a spectrum
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectraMerger &>())
        .def("__copy__", [](const OpenMS::SpectraMerger& self) { return OpenMS::SpectraMerger(self); })
        .def("__deepcopy__", [](const OpenMS::SpectraMerger& self, nb::dict) { return OpenMS::SpectraMerger(self); }, "memo"_a)
        .def("mergeSpectraBlockWise", [](OpenMS::SpectraMerger& self, OpenMS::MSExperiment& exp) { self.mergeSpectraBlockWise(exp); }, "exp"_a, "Merges spectra block-wise")
        .def("mergeSpectraPrecursors", [](OpenMS::SpectraMerger& self, OpenMS::MSExperiment& exp) { self.mergeSpectraPrecursors(exp); }, "exp"_a, "Merges spectra with similar precursors")
        .def("average", [](OpenMS::SpectraMerger& self, OpenMS::MSExperiment& exp, const OpenMS::String& average_type, int ms_level) { self.average(exp, average_type, ms_level); }, "exp"_a, "average_type"_a, "ms_level"_a = -1, "Averages spectra")
        ;

    // -----------------------------------------------------------------------
    // SpectralDeconvolution
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectralDeconvolution, OpenMS::DefaultParamHandler>(m, "SpectralDeconvolution", 
        R"doc(
DefaultParamHandler

Spectral deconvolution algorithm for top-down MS.
From MSSpectrum, this class outputs DeconvolvedSpectrum.
Deconvolution takes three steps:
1) decharging and select candidate masses - speed up via binning
2) collecting isotopes from the candidate masses and deisotoping
3) scoring and filter out low scoring masses
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectralDeconvolution &>())
        .def("__copy__", [](const OpenMS::SpectralDeconvolution& self) { return OpenMS::SpectralDeconvolution(self); })
        .def("__deepcopy__", [](const OpenMS::SpectralDeconvolution& self, nb::dict) { return OpenMS::SpectralDeconvolution(self); }, "memo"_a)
        .def("performSpectrumDeconvolution", [](OpenMS::SpectralDeconvolution& self, const OpenMS::MSSpectrum& spec, int scan_number, const OpenMS::PeakGroup& precursor_peak_group) { return self.performSpectrumDeconvolution(spec, scan_number, precursor_peak_group); }, "spec"_a, "scan_number"_a, "precursor_peak_group"_a)
        .def("getDeconvolvedSpectrum", [](OpenMS::SpectralDeconvolution& self) -> OpenMS::DeconvolvedSpectrum & { return self.getDeconvolvedSpectrum(); }, nb::rv_policy::reference_internal, 
            R"doc(
Main deconvolution function that generates the deconvolved spectrum.
:param spec: The original spectrum
:param scan_number: Scan number from input spectrum
:param precursor_peak_group: Precursor peak group (for MS2+)
Result access
)doc")
        .def("getAveragine", [](OpenMS::SpectralDeconvolution& self) -> const OpenMS::FLASHHelperClasses::PrecalculatedAveragine & { return self.getAveragine(); }, nb::rv_policy::reference_internal, "Return the deconvolved spectrum after performSpectrumDeconvolution is called")
        .def("setAveragine", [](OpenMS::SpectralDeconvolution& self, const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& avg) { return self.setAveragine(avg); }, "avg"_a, "Get calculated averagine. Call after calculateAveragine is called.")
        .def("setTargetMasses", [](OpenMS::SpectralDeconvolution& self, const std::vector<double>& masses, bool exclude) { return self.setTargetMasses(masses, exclude); }, "masses"_a, "exclude"_a = false, "Set the precalculated averagine")
        .def("calculateAveragine", [](OpenMS::SpectralDeconvolution& self, bool use_RNA_averagine) { return self.calculateAveragine(use_RNA_averagine); }, "use_RNA_averagine"_a, "Set targeted or excluded masses for targeted deconvolution. Masses are targeted or excluded in all ms levels.")
        .def("setToleranceEstimation", [](OpenMS::SpectralDeconvolution& self) { return self.setToleranceEstimation(); }, "Precalculate averagine (for predefined mass bins) to speed up averagine generation")
        .def_static("getNominalMass", [](double mass) { return OpenMS::SpectralDeconvolution::getNominalMass(mass); }, "mass"_a, "Set target decoy type for the SpectralDeconvolution run")
        .def_static("getCosine", [](const std::vector<float>& a, int a_start, int a_end, const OpenMS::IsotopeDistribution& b, int offset, int min_iso_len) { return OpenMS::SpectralDeconvolution::getCosine(a, a_start, a_end, b, offset, min_iso_len); }, "a"_a, "a_start"_a, "a_end"_a, "b"_a, "offset"_a, "min_iso_len"_a, "Convert double mass to nominal mass (integer)")
        .def_static("getIsotopeCosineAndIsoOffset", [](double mono_mass, const std::vector<float>& per_isotope_intensities, const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& avg, int iso_int_shift, int window_width, const std::vector<double>& excluded_masses) {
            int offset = 0;
            auto cosine = OpenMS::SpectralDeconvolution::getIsotopeCosineAndIsoOffset(mono_mass, per_isotope_intensities, offset, avg, iso_int_shift, window_width, excluded_masses);
            return nb::make_tuple(cosine, offset);
        }, "mono_mass"_a, "per_isotope_intensities"_a, "avg"_a, "iso_int_shift"_a, "window_width"_a, "excluded_masses"_a, "Calculate cosine between two vectors with optimization parameters")
        .def("setTargetDecoyType", [](OpenMS::SpectralDeconvolution& self, OpenMS::PeakGroup::TargetDecoyType target_decoy_type, const OpenMS::DeconvolvedSpectrum& target_dspec_for_decoy_calcualtion) { return self.setTargetDecoyType(target_decoy_type, target_dspec_for_decoy_calcualtion); }, "target_decoy_type"_a, "target_dspec_for_decoy_calcualtion"_a, "When estimating tolerance, set max_mass_dalton_tolerance to a large value")
        ;

    // -----------------------------------------------------------------------
    // Spectrum
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Interfaces::Spectrum>(m, "Spectrum", "_Interfaces_BinaryDataArray := BinaryDataArray")
        .def(nb::init<>())
        .def_rw("defaultArrayLength", &OpenMS::Interfaces::Spectrum::defaultArrayLength)

        .def("setMZArray", [](OpenMS::Interfaces::Spectrum& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setMZArray(arr);
        }, "data"_a, "Set m/z array from list")

        .def("setIntensityArray", [](OpenMS::Interfaces::Spectrum& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")

        .def("getMZArray", [](const OpenMS::Interfaces::Spectrum& self) {
            auto arr = self.getMZArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get m/z array")

        .def("getIntensityArray", [](const OpenMS::Interfaces::Spectrum& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array")
        ;

    // -----------------------------------------------------------------------
    // SpectrumAlignment
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumAlignment, OpenMS::DefaultParamHandler>(m, "SpectrumAlignment", 
        R"doc(
DefaultParamHandler

Aligns the peaks of two sorted spectra.
Method 1: Using a banded (width via 'tolerance' parameter) alignment if absolute tolerances are given.
Scoring function is the m/z distance between peaks. Intensity does not play a role!
Method 2: If relative tolerance (ppm) is specified a simple matching of peaks is performed:
Peaks from s1 (usually the theoretical spectrum) are assigned to the closest peak in s2 if it lies in the tolerance window.
Note: A peak in s2 can be matched to none, one or multiple peaks in s1. Peaks in s1 may be matched to none or one peak in s2.
Note: Intensity is ignored.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectrumAlignment &>())
        .def("__copy__", [](const OpenMS::SpectrumAlignment& self) { return OpenMS::SpectrumAlignment(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumAlignment& self, nb::dict) { return OpenMS::SpectrumAlignment(self); }, "memo"_a)
        .def("getSpectrumAlignment", [](const OpenMS::SpectrumAlignment& self, const OpenMS::MSSpectrum& s1, const OpenMS::MSSpectrum& s2) {
            std::vector<std::pair<OpenMS::Size, OpenMS::Size>> alignment;
            self.getSpectrumAlignment(alignment, s1, s2);
            nb::list result;
            for (const auto& p : alignment) {
                result.append(nb::make_tuple(p.first, p.second));
            }
            return result;
        }, "s1"_a, "s2"_a, "Aligns two spectra and returns a list of (index_s1, index_s2) tuples")
        ;

    // -----------------------------------------------------------------------
    // SpectrumAnnotator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumAnnotator, OpenMS::DefaultParamHandler>(m, "SpectrumAnnotator", 
        R"doc(
Annotates spectra from identifications and theoretical spectra or
identifications from spectra and theoretical spectra matching with
various options
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectrumAnnotator &>())
        .def("__copy__", [](const OpenMS::SpectrumAnnotator& self) { return OpenMS::SpectrumAnnotator(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumAnnotator& self, nb::dict) { return OpenMS::SpectrumAnnotator(self); }, "memo"_a)
        .def("annotateMatches", [](const OpenMS::SpectrumAnnotator& self, OpenMS::MSSpectrum& spec, const OpenMS::PeptideHit& ph, const OpenMS::TheoreticalSpectrumGenerator& tg, const OpenMS::SpectrumAlignment& sa) { return self.annotateMatches(spec, ph, tg, sa); }, "spec"_a, "ph"_a, "tg"_a, "sa"_a)
        .def("addIonMatchStatistics", [](const OpenMS::SpectrumAnnotator& self, OpenMS::PeptideIdentification& pi, OpenMS::MSSpectrum& spec, const OpenMS::TheoreticalSpectrumGenerator& tg, const OpenMS::SpectrumAlignment& sa) { return self.addIonMatchStatistics(pi, spec, tg, sa); }, "pi"_a, "spec"_a, "tg"_a, "sa"_a, 
            R"doc(
Adds ion match annotation to the `spec` input spectrum
:param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
:param ph: A spectrum identifications to be used for the annotation, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
:param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
:param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
)doc")
        .def("addPeakAnnotationsToPeptideHit", [](const OpenMS::SpectrumAnnotator& self, OpenMS::PeptideHit& ph, const OpenMS::MSSpectrum& spec, const OpenMS::TheoreticalSpectrumGenerator& tg, const OpenMS::SpectrumAlignment& sa, bool include_unmatched_peaks) { return self.addPeakAnnotationsToPeptideHit(ph, spec, tg, sa, include_unmatched_peaks); }, "ph"_a, "spec"_a, "tg"_a, "sa"_a, "include_unmatched_peaks"_a = false, 
            R"doc(
Adds ion match statistics to `pi` PeptideIdentifcation
:param pi: A spectrum identifications to be annotated, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
:param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
:param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
:param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
)doc")
        ;

    // -----------------------------------------------------------------------
    // SqrtScaler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SqrtScaler, OpenMS::DefaultParamHandler>(m, "SqrtScaler", 
        R"doc(
DefaultParamHandler

Scales the intensity of peaks to the sqrt
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SqrtScaler &>())
        .def("__copy__", [](const OpenMS::SqrtScaler& self) { return OpenMS::SqrtScaler(self); })
        .def("__deepcopy__", [](const OpenMS::SqrtScaler& self, nb::dict) { return OpenMS::SqrtScaler(self); }, "memo"_a)
        .def("filterPeakSpectrum", [](OpenMS::SqrtScaler& self, OpenMS::MSSpectrum& spectrum) { return self.filterPeakSpectrum(spectrum); }, "spectrum"_a)
        .def("filterPeakMap", [](OpenMS::SqrtScaler& self, OpenMS::MSExperiment& exp) { return self.filterPeakMap(exp); }, "exp"_a)
        ;

    // -----------------------------------------------------------------------
    // StablePairFinder
    // -----------------------------------------------------------------------
    auto stablepairfinder_class = nb::class_<OpenMS::StablePairFinder, OpenMS::BaseGroupFinder>(m, "StablePairFinder", 
        R"doc(
This class implements a pair finding algorithm for consensus features
BaseGroupFinder
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::StablePairFinder& self, const std::vector<OpenMS::ConsensusMap>& input_maps, OpenMS::ConsensusMap& result_map) { self.run(input_maps, result_map); }, "input_maps"_a, "result_map"_a)
        ;
    def_ProgressLogger<OpenMS::StablePairFinder>(stablepairfinder_class);

    // -----------------------------------------------------------------------
    // SwathFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SwathFile, OpenMS::ProgressLogger>(m, "SwathFile", 
        R"doc(
File adapter for Swath files. * * This class can load SWATH files in
different storage versions. The most * convenient file is a single
MzML file which contains one experiment. * However, also the loading
of a list of files is supported (loadSplit) * where it is assumed that
each individual file only contains scans from one * precursor
isolation window (one SWATH). Finally, experimental support for *
mzXML is available but needs to be selected with a specific compile
flag * (this is not for everyday use). *
ProgressLogger
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SwathFile &>())
        .def("__copy__", [](const OpenMS::SwathFile& self) { return OpenMS::SwathFile(self); })
        .def("__deepcopy__", [](const OpenMS::SwathFile& self, nb::dict) { return OpenMS::SwathFile(self); }, "memo"_a)
        .def("loadSplit", [](OpenMS::SwathFile& self, std::vector<OpenMS::String> file_list, const OpenMS::String& tmp, std::shared_ptr<OpenMS::ExperimentalSettings>& exp_meta, const OpenMS::String& readoptions) { return self.loadSplit(file_list, tmp, exp_meta, readoptions); }, "file_list"_a, "tmp"_a, "exp_meta"_a, "readoptions"_a = "normal", "Loads a Swath run from a list of split mzML files")
        .def("loadMzML", [](OpenMS::SwathFile& self, const OpenMS::String& file, const OpenMS::String& tmp, std::shared_ptr<OpenMS::ExperimentalSettings>& exp_meta, const OpenMS::String& readoptions, OpenMS::Interfaces::IMSDataConsumer * plugin_consumer) { return self.loadMzML(file, tmp, exp_meta, readoptions, plugin_consumer); }, "file"_a, "tmp"_a, "exp_meta"_a, "readoptions"_a = "normal", "plugin_consumer"_a)
        .def("loadMzXML", [](OpenMS::SwathFile& self, const OpenMS::String& file, const OpenMS::String& tmp, std::shared_ptr<OpenMS::ExperimentalSettings>& exp_meta, const OpenMS::String& readoptions) { return self.loadMzXML(file, tmp, exp_meta, readoptions); }, "file"_a, "tmp"_a, "exp_meta"_a, "readoptions"_a = "normal", "Loads a Swath run from a single mzXML file")
        ;

    // -----------------------------------------------------------------------
    // SwathMapMassCorrection
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SwathMapMassCorrection, OpenMS::DefaultParamHandler>(m, "SwathMapMassCorrection", 
        R"doc(
A class containing correction functions for Swath MS maps * * This
class can use a set of pre-determined points in a Swath-MS map to *
correct all maps according to the m/z shift found in those fixed
points. *
)doc")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // TMTEighteenPlexQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TMTEighteenPlexQuantitationMethod, OpenMS::IsobaricQuantitationMethod>(m, "TMTEighteenPlexQuantitationMethod", 
        R"doc(
TMT 18-plex isobaric labeling quantitation method. Provides channel information
and isotope correction matrix for Thermo Scientific TMTpro 18-plex reagents
which extends TMTpro 16-plex with two additional channels
IsobaricQuantitationMethod
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TMTEighteenPlexQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::TMTEighteenPlexQuantitationMethod& self) { return OpenMS::TMTEighteenPlexQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::TMTEighteenPlexQuantitationMethod& self, nb::dict) { return OpenMS::TMTEighteenPlexQuantitationMethod(self); }, "memo"_a)
        .def("getChannelInformation", [](const OpenMS::TMTEighteenPlexQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::TMTEighteenPlexQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::TMTEighteenPlexQuantitationMethod& self) { return self.getIsotopeCorrectionMatrix(); }, "Returns the isotope correction matrix for correcting reporter ion intensities")
        .def("getReferenceChannel", [](const OpenMS::TMTEighteenPlexQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")
        ;

    // -----------------------------------------------------------------------
    // TMTElevenPlexQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TMTElevenPlexQuantitationMethod, OpenMS::IsobaricQuantitationMethod>(m, "TMTElevenPlexQuantitationMethod", 
        R"doc(
TMT 11-plex isobaric labeling quantitation method. Provides channel information
and isotope correction matrix for Thermo Scientific TMT 11-plex reagents
which extends TMT 10-plex with an additional 131C channel
IsobaricQuantitationMethod
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TMTElevenPlexQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::TMTElevenPlexQuantitationMethod& self) { return OpenMS::TMTElevenPlexQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::TMTElevenPlexQuantitationMethod& self, nb::dict) { return OpenMS::TMTElevenPlexQuantitationMethod(self); }, "memo"_a)
        .def("getChannelInformation", [](const OpenMS::TMTElevenPlexQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::TMTElevenPlexQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::TMTElevenPlexQuantitationMethod& self) { return self.getIsotopeCorrectionMatrix(); }, "Returns the isotope correction matrix for correcting reporter ion intensities")
        .def("getReferenceChannel", [](const OpenMS::TMTElevenPlexQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")
        ;

    // -----------------------------------------------------------------------
    // TMTSixPlexQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TMTSixPlexQuantitationMethod, OpenMS::IsobaricQuantitationMethod>(m, "TMTSixPlexQuantitationMethod", 
        R"doc(
TMT 6-plex isobaric labeling quantitation method. Provides channel information
and isotope correction matrix for Thermo Scientific TMT 6-plex reagents
with reporter ions at 126, 127, 128, 129, 130, and 131 m/z
IsobaricQuantitationMethod
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TMTSixPlexQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::TMTSixPlexQuantitationMethod& self) { return OpenMS::TMTSixPlexQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::TMTSixPlexQuantitationMethod& self, nb::dict) { return OpenMS::TMTSixPlexQuantitationMethod(self); }, "memo"_a)
        .def("getChannelInformation", [](const OpenMS::TMTSixPlexQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::TMTSixPlexQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::TMTSixPlexQuantitationMethod& self) { return self.getIsotopeCorrectionMatrix(); }, "Returns the isotope correction matrix for correcting reporter ion intensities")
        .def("getReferenceChannel", [](const OpenMS::TMTSixPlexQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")
        ;

    // -----------------------------------------------------------------------
    // TMTSixteenPlexQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TMTSixteenPlexQuantitationMethod, OpenMS::IsobaricQuantitationMethod>(m, "TMTSixteenPlexQuantitationMethod", 
        R"doc(
TMT 16-plex isobaric labeling quantitation method. Provides channel information
and isotope correction matrix for Thermo Scientific TMTpro 16-plex reagents
with reporter ions spanning the 126-134 m/z range
IsobaricQuantitationMethod
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TMTSixteenPlexQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::TMTSixteenPlexQuantitationMethod& self) { return OpenMS::TMTSixteenPlexQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::TMTSixteenPlexQuantitationMethod& self, nb::dict) { return OpenMS::TMTSixteenPlexQuantitationMethod(self); }, "memo"_a)
        .def("getChannelInformation", [](const OpenMS::TMTSixteenPlexQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::TMTSixteenPlexQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::TMTSixteenPlexQuantitationMethod& self) { return self.getIsotopeCorrectionMatrix(); }, "Returns the isotope correction matrix for correcting reporter ion intensities")
        .def("getReferenceChannel", [](const OpenMS::TMTSixteenPlexQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")
        ;

    // -----------------------------------------------------------------------
    // TMTTenPlexQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TMTTenPlexQuantitationMethod, OpenMS::IsobaricQuantitationMethod>(m, "TMTTenPlexQuantitationMethod", 
        R"doc(
TMT 10-plex isobaric labeling quantitation method. Provides channel information
and isotope correction matrix for Thermo Scientific TMT 10-plex reagents
with reporter ions spanning the 126-131 m/z range with additional mass variants
IsobaricQuantitationMethod
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TMTTenPlexQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::TMTTenPlexQuantitationMethod& self) { return OpenMS::TMTTenPlexQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::TMTTenPlexQuantitationMethod& self, nb::dict) { return OpenMS::TMTTenPlexQuantitationMethod(self); }, "memo"_a)
        .def("getChannelInformation", [](const OpenMS::TMTTenPlexQuantitationMethod& self) -> const std::vector<OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation> & { return self.getChannelInformation(); }, nb::rv_policy::reference_internal, "Returns information on the different channels used by this quantitation method")
        .def("getNumberOfChannels", [](const OpenMS::TMTTenPlexQuantitationMethod& self) { return self.getNumberOfChannels(); }, "Returns the number of channels available for this quantitation method")
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::TMTTenPlexQuantitationMethod& self) { return self.getIsotopeCorrectionMatrix(); }, "Returns the isotope correction matrix for correcting reporter ion intensities")
        .def("getReferenceChannel", [](const OpenMS::TMTTenPlexQuantitationMethod& self) { return self.getReferenceChannel(); }, "Returns the index of the reference channel used for ratio calculation")
        ;

    // -----------------------------------------------------------------------
    // TSE_Match (TargetedSpectraExtractor::Match)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedSpectraExtractor::Match>(m, "TSE_Match",
        "A match between a spectrum and a score from TargetedSpectraExtractor")
        .def(nb::init<>())
        .def(nb::init<OpenMS::MSSpectrum, double>(), "spectrum"_a, "score"_a)
        .def("__copy__", [](const OpenMS::TargetedSpectraExtractor::Match& self) { return OpenMS::TargetedSpectraExtractor::Match(self); })
        .def("__deepcopy__", [](const OpenMS::TargetedSpectraExtractor::Match& self, nb::dict) { return OpenMS::TargetedSpectraExtractor::Match(self); }, "memo"_a)
        .def_rw("spectrum", &OpenMS::TargetedSpectraExtractor::Match::spectrum)
        .def_rw("score", &OpenMS::TargetedSpectraExtractor::Match::score)
        ;

    // -----------------------------------------------------------------------
    // TargetedSpectraExtractor
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedSpectraExtractor, OpenMS::DefaultParamHandler>(m, "TargetedSpectraExtractor", 
        R"doc(
DefaultParamHandler

Filter, annotate, pick, and score spectra based on a target list
This class processes spectra from DDA experiments against a target transition list.
It provides a complete pipeline from raw spectra to scored, selected spectra.
Workflow:
1. annotateSpectra() - Match spectra to targets by precursor MZ and RT
2. pickSpectrum() - Smooth and pick peaks
3. scoreSpectra() - Score by TIC, FWHM, SNR with configurable weights
4. selectSpectra() - Choose the best spectrum per transition
Or use extractSpectra() for the full pipeline in one call.
Example usage:
.. code-block:: python
from pyopenms import *
# Load experiment and targets
exp = MSExperiment()
MzMLFile().load("experiment.mzML", exp)
targets = TargetedExperiment()
TraMLFile().load("targets.traML", targets)
# Extract matching spectra
extractor = TargetedSpectraExtractor()
extracted, features = extractor.extractSpectra(exp, targets, True)
print(f"Found {len(extracted)} matching spectra")
)doc")
        .def(nb::init<>())
        .def("annotateSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const std::vector<OpenMS::MSSpectrum>& spectra, const OpenMS::TargetedExperiment& targeted_exp, OpenMS::FeatureMap& features, bool compute_features) {
            std::vector<OpenMS::MSSpectrum> annotated_spectra;
            self.annotateSpectra(spectra, targeted_exp, annotated_spectra, features, compute_features);
            return annotated_spectra;
        }, "spectra"_a, "targeted_exp"_a, "features"_a, "compute_features"_a, "Get default parameters for the extractor")
        .def("annotateSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const std::vector<OpenMS::MSSpectrum>& spectra, const OpenMS::TargetedExperiment& targeted_exp) {
            std::vector<OpenMS::MSSpectrum> annotated_spectra;
            self.annotateSpectra(spectra, targeted_exp, annotated_spectra);
            return annotated_spectra;
        }, "spectra"_a, "targeted_exp"_a, 
            R"doc(
Annotate spectra that match target transitions
Filters spectra that fall within the precursor RT window and MZ tolerance.
:param spectra: Input spectra to filter
:param targeted_exp: Target transition list
:param annotated_spectra: Output annotated spectra
:param features: Output FeatureMap with transition info
:param compute_features: If True, populate features output
)doc")
        .def("annotateSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const std::vector<OpenMS::MSSpectrum>& spectra, const OpenMS::FeatureMap& ms1_features, OpenMS::FeatureMap& ms2_features) {
            std::vector<OpenMS::MSSpectrum> annotated_spectra;
            self.annotateSpectra(spectra, ms1_features, ms2_features, annotated_spectra);
            return annotated_spectra;
        }, "spectra"_a, "ms1_features"_a, "ms2_features"_a, 
            R"doc(
Annotate spectra that match target transitions (without features)
:param spectra: Input spectra to filter
:param targeted_exp: Target transition list
:param annotated_spectra: Output annotated spectra
)doc")
        .def("searchSpectrum", [](const OpenMS::TargetedSpectraExtractor& self, OpenMS::FeatureMap& feat_map, OpenMS::FeatureMap& feat_map_output, bool add_unidentified_features) { return self.searchSpectrum(feat_map, feat_map_output, add_unidentified_features); }, "feat_map"_a, "feat_map_output"_a, "add_unidentified_features"_a = false, 
            R"doc(
Annotate spectra using MS1 and MS2 feature maps
:param spectra: Input spectra
:param ms1_features: MS1 feature map
:param ms2_features: MS2 feature map
:param annotated_spectra: Output annotated spectra
)doc")
        .def("pickSpectrum", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::MSSpectrum& spectrum, OpenMS::MSSpectrum& picked_spectrum) { return self.pickSpectrum(spectrum, picked_spectrum); }, "spectrum"_a, "picked_spectrum"_a, 
            R"doc(
Search for matching features between MS1 and MS2
:param ms1_features: MS1 feature map
:param ms2_features: MS2 feature map (modified in place)
:param add_unidentified_features: Add features without matches
)doc")
        .def("scoreSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const std::vector<OpenMS::MSSpectrum>& annotated_spectra, const std::vector<OpenMS::MSSpectrum>& picked_spectra, OpenMS::FeatureMap& features, bool compute_features) {
            std::vector<OpenMS::MSSpectrum> scored_spectra;
            self.scoreSpectra(annotated_spectra, picked_spectra, features, scored_spectra, compute_features);
            return scored_spectra;
        }, "annotated_spectra"_a, "picked_spectra"_a, "features"_a, "compute_features"_a, 
            R"doc(
Smooth and pick peaks from a spectrum
Applies Gaussian or Savitzky-Golay smoothing followed by peak picking.
:param spectrum: Input spectrum
:param picked_spectrum: Output picked spectrum
)doc")
        .def("scoreSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const std::vector<OpenMS::MSSpectrum>& annotated_spectra, const std::vector<OpenMS::MSSpectrum>& picked_spectra) {
            std::vector<OpenMS::MSSpectrum> scored_spectra;
            self.scoreSpectra(annotated_spectra, picked_spectra, scored_spectra);
            return scored_spectra;
        }, "annotated_spectra"_a, "picked_spectra"_a, 
            R"doc(
Score spectra based on TIC, FWHM, and SNR
Scores are computed using configurable weights for each metric.
:param annotated_spectra: Annotated input spectra
:param picked_spectra: Peak-picked spectra
:param features: Feature map with transition info
:param scored_spectra: Output scored spectra
)doc")
        .def("selectSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const std::vector<OpenMS::MSSpectrum>& scored_spectra, const OpenMS::FeatureMap& features, OpenMS::FeatureMap& selected_features, bool compute_features) {
            std::vector<OpenMS::MSSpectrum> selected_spectra;
            self.selectSpectra(scored_spectra, features, selected_spectra, selected_features, compute_features);
            return selected_spectra;
        }, "scored_spectra"_a, "features"_a, "selected_features"_a, "compute_features"_a, 
            R"doc(
Score spectra (without feature map)
:param annotated_spectra: Annotated input spectra
:param picked_spectra: Peak-picked spectra
:param scored_spectra: Output scored spectra
)doc")
        .def("selectSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const std::vector<OpenMS::MSSpectrum>& scored_spectra) {
            std::vector<OpenMS::MSSpectrum> selected_spectra;
            self.selectSpectra(scored_spectra, selected_spectra);
            return selected_spectra;
        }, "scored_spectra"_a, 
            R"doc(
Select the best spectrum for each transition
:param scored_spectra: Scored input spectra
:param features: Feature map with transition info
:param selected_spectra: Output selected spectra
:param selected_features: Output selected features
)doc")
        .def("extractSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::MSExperiment& experiment, const OpenMS::TargetedExperiment& targeted_exp, bool compute_features) { std::vector<OpenMS::MSSpectrum> extracted_spectra; OpenMS::FeatureMap extracted_features; self.extractSpectra(experiment, targeted_exp, extracted_spectra, extracted_features, compute_features); return nb::make_tuple(extracted_spectra, extracted_features); }, "experiment"_a, "targeted_exp"_a, "compute_features"_a, 
            R"doc(
Select the best spectrum for each transition (without features)
:param scored_spectra: Scored input spectra
:param selected_spectra: Output selected spectra
)doc")
        .def("extractSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::MSExperiment& experiment, const OpenMS::TargetedExperiment& targeted_exp) {
            std::vector<OpenMS::MSSpectrum> extracted_spectra;
            self.extractSpectra(experiment, targeted_exp, extracted_spectra);
            return extracted_spectra;
        }, "experiment"_a, "targeted_exp"_a, 
            R"doc(
Complete pipeline: annotate, pick, score, and select spectra
Runs the full extraction workflow in a single call.
:param experiment: Input MS experiment
:param targeted_exp: Target transition list
:param extracted_spectra: Output extracted spectra
:param features: Output feature map
:param compute_features: If True, populate features
)doc")
        .def("extractSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::MSExperiment& experiment, const OpenMS::FeatureMap& ms1_features) {
            std::vector<OpenMS::MSSpectrum> extracted_spectra;
            self.extractSpectra(experiment, ms1_features, extracted_spectra);
            return extracted_spectra;
        }, "experiment"_a, "ms1_features"_a, 
            R"doc(
Complete pipeline: annotate, pick, score, and select spectra
Runs the full extraction workflow in a single call.
:param experiment: Input MS experiment
:param targeted_exp: Target transition list
:param extracted_spectra: Output extracted spectra
:param features: Output feature map
:param compute_features: If True, populate features
)doc")
        .def("extractSpectra", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::MSExperiment& experiment, const OpenMS::FeatureMap& ms1_features, OpenMS::FeatureMap& extracted_features) {
            std::vector<OpenMS::MSSpectrum> extracted_spectra;
            self.extractSpectra(experiment, ms1_features, extracted_spectra, extracted_features);
            return extracted_spectra;
        }, "experiment"_a, "ms1_features"_a, "extracted_features"_a, 
            R"doc(
Select the best spectrum for each transition (without features)
:param scored_spectra: Scored input spectra
:param selected_spectra: Output selected spectra
)doc")
        .def("constructTransitionsList", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::FeatureMap& ms1_features, const OpenMS::FeatureMap& ms2_features, OpenMS::TargetedExperiment& t_exp) { return self.constructTransitionsList(ms1_features, ms2_features, t_exp); }, "ms1_features"_a, "ms2_features"_a, "t_exp"_a, 
            R"doc(
Extract spectra using MS2 features
:param experiment: Input MS experiment
:param ms2_features: MS2 feature map
:param extracted_spectra: Output extracted spectra
)doc")
        .def("storeSpectraMSP", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::String& filename, OpenMS::MSExperiment& experiment) { return self.storeSpectraMSP(filename, experiment); }, "filename"_a, "experiment"_a, 
            R"doc(
Construct a transition list from MS1 and MS2 features
:param ms1_features: MS1 feature map
:param ms2_features: MS2 feature map
:param t_exp: Output TargetedExperiment
)doc")
        .def("mergeFeatures", [](const OpenMS::TargetedSpectraExtractor& self, const OpenMS::FeatureMap& fmap_input, OpenMS::FeatureMap& fmap_output) { return self.mergeFeatures(fmap_input, fmap_output); }, "fmap_input"_a, "fmap_output"_a, 
            R"doc(
Store spectra in MSP format
:param filename: Output filename
:param experiment: Experiment to store
)doc")

        .def("getDefaultParameters", [](const OpenMS::TargetedSpectraExtractor& self, OpenMS::Param& params) {
            self.getDefaultParameters(params);
        }, "params"_a, "Get default parameters (fills output param)")
        .def("getDefaultParameters", [](const OpenMS::TargetedSpectraExtractor& self) {
            OpenMS::Param params;
            self.getDefaultParameters(params);
            return params;
        }, "Get default parameters (returns Param)")
        ;

    // -----------------------------------------------------------------------
    // TextFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TextFile>(m, "TextFile", "OpenMS class TextFile")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::TextFile& self) { return OpenMS::TextFile(self); })
        .def("__deepcopy__", [](const OpenMS::TextFile& self, nb::dict) { return OpenMS::TextFile(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, bool, int, bool, OpenMS::String>())
        .def("load", [](OpenMS::TextFile& self, const OpenMS::String& filename, bool trim_lines, int first_n, bool skip_empty_lines, const OpenMS::String& comment_symbol) { return self.load(filename, trim_lines, first_n, skip_empty_lines, comment_symbol); }, "filename"_a, "trim_lines"_a = false, "first_n"_a = -1, "skip_empty_lines"_a = false, "comment_symbol"_a = "")
        .def("store", [](OpenMS::TextFile& self, const OpenMS::String& filename) { return self.store(filename); }, "filename"_a, "Writes the data to a file")
        .def("addLine", [](OpenMS::TextFile& self, const OpenMS::String& line) { self.addLine(line); }, "line"_a, "Appends a line to the internal buffer")
        ;

    // -----------------------------------------------------------------------
    // CsvFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CsvFile, OpenMS::TextFile>(m, "CsvFile", 
        R"doc(
This class handles csv files. Currently only loading is implemented.
Does NOT support comment lines!
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::CsvFile& self) { return OpenMS::CsvFile(self); })
        .def("__deepcopy__", [](const OpenMS::CsvFile& self, nb::dict) { return OpenMS::CsvFile(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, char, bool, int>())
        .def("load", [](OpenMS::CsvFile& self, const OpenMS::String& filename, char is, bool ie, int first_n) { return self.load(filename, is, ie, first_n); }, "filename"_a, "is"_a = ',', "ie"_a = false, "first_n"_a = -1, "Loads data from a text file")
        .def("store", [](OpenMS::CsvFile& self, const OpenMS::String& filename) { return self.store(filename); }, "filename"_a, "Stores the buffer's content into a file")
        .def("addRow", [](OpenMS::CsvFile& self, const std::vector<OpenMS::String>& list) { return self.addRow(list); }, "list"_a, "Add a row to the buffer")
        .def("clear", [](OpenMS::CsvFile& self) { return self.clear(); }, "Clears the buffer")
        .def("getRow", [](const OpenMS::CsvFile& self, size_t row) {
            std::vector<OpenMS::String> list;
            auto result = self.getRow(row, list);
            return nb::make_tuple(result, list);
        }, "row"_a, "Writes all items from a row to list")
        ;

    // -----------------------------------------------------------------------
    // AbsoluteQuantitationMethodFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitationMethodFile, OpenMS::CsvFile>(m, "AbsoluteQuantitationMethodFile", 
        R"doc(
File adapter for AbsoluteQuantitationMethod files. Loads and stores .csv or .tsv
files describing absolute quantitation methods including calibration curve
parameters, limits of detection/quantitation, and transformation models
)doc")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // MRMFeaturePickerFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeaturePickerFile, OpenMS::CsvFile>(m, "MRMFeaturePickerFile", 
        R"doc(
_MRMFeaturePickerFile_ loads components and components groups parameters from a .csv file
The structures defined in [MRMFeaturePicker](@ref MRMFeaturePicker) are used
It is required that columns `component_name` and `component_group_name` are present.
Lines whose `component_name`'s or `component_group_name`'s value is an empty string, will be skipped.
The class supports the absence of information within other columns.
A reduced example of the expected format (fewer columns are shown here):
> component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length
> arg-L.arg-L_1.Heavy,arg-L,2,15
> arg-L.arg-L_1.Light,arg-L,2,17
> orn.orn_1.Heavy,orn,3,21
> orn.orn_1.Light,orn,3,13
)doc")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // MRMFeatureQCFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureQCFile, OpenMS::CsvFile>(m, "MRMFeatureQCFile", 
        R"doc(
File adapter for MRMFeatureQC files
Loads and stores .csv or .tsv files describing an MRMFeatureQC
)doc")
        .def(nb::init<>())
        .def("load", [](const OpenMS::MRMFeatureQCFile& self, const OpenMS::String& filename, bool is_component_group) { OpenMS::MRMFeatureQC mrmfqc; self.load(filename, mrmfqc, is_component_group); return mrmfqc; }, "filename"_a, "is_component_group"_a)
        ;

    // -----------------------------------------------------------------------
    // TheoreticalSpectrumGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TheoreticalSpectrumGenerator, OpenMS::DefaultParamHandler>(m, "TheoreticalSpectrumGenerator", 
        R"doc(
Generates theoretical spectra for peptides with various options
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TheoreticalSpectrumGenerator &>())
        .def("__copy__", [](const OpenMS::TheoreticalSpectrumGenerator& self) { return OpenMS::TheoreticalSpectrumGenerator(self); })
        .def("__deepcopy__", [](const OpenMS::TheoreticalSpectrumGenerator& self, nb::dict) { return OpenMS::TheoreticalSpectrumGenerator(self); }, "memo"_a)

        .def("getSpectrum", [](OpenMS::TheoreticalSpectrumGenerator& self, OpenMS::MSSpectrum& spec, const OpenMS::AASequence& peptide, int min_charge, int max_charge) {
            self.getSpectrum(spec, peptide, min_charge, max_charge);
        }, "spec"_a, "peptide"_a, "min_charge"_a, "max_charge"_a, "Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters")
        ;

    // -----------------------------------------------------------------------
    // TheoreticalSpectrumGeneratorXLMS
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TheoreticalSpectrumGeneratorXLMS, OpenMS::DefaultParamHandler>(m, "TheoreticalSpectrumGeneratorXLMS", 
        R"doc(
Generates theoretical spectra for cross-linked peptides
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TheoreticalSpectrumGeneratorXLMS &>())
        .def("__copy__", [](const OpenMS::TheoreticalSpectrumGeneratorXLMS& self) { return OpenMS::TheoreticalSpectrumGeneratorXLMS(self); })
        .def("__deepcopy__", [](const OpenMS::TheoreticalSpectrumGeneratorXLMS& self, nb::dict) { return OpenMS::TheoreticalSpectrumGeneratorXLMS(self); }, "memo"_a)
        .def("getLinearIonSpectrum", [](const OpenMS::TheoreticalSpectrumGeneratorXLMS& self, size_t link_pos, bool frag_alpha, int charge, size_t link_pos_2) { OpenMS::MSSpectrum spectrum; OpenMS::AASequence peptide; self.getLinearIonSpectrum(spectrum, peptide, link_pos, frag_alpha, charge, link_pos_2); return std::make_tuple(spectrum, peptide); }, "link_pos"_a, "frag_alpha"_a, "charge"_a, "link_pos_2"_a)
        .def("getXLinkIonSpectrum", [](const OpenMS::TheoreticalSpectrumGeneratorXLMS& self, bool frag_alpha, int mincharge, int maxcharge) { OpenMS::MSSpectrum spectrum; OpenMS::OPXLDataStructs::ProteinProteinCrossLink crosslink; self.getXLinkIonSpectrum(spectrum, crosslink, frag_alpha, mincharge, maxcharge); return std::make_tuple(spectrum, crosslink); }, "frag_alpha"_a, "mincharge"_a, "maxcharge"_a, 
            R"doc(
Generates fragment ions containing the cross-linker for one peptide
B-ions are generated from the first linked position up to the end of the peptide,
y-ions are generated from the beginning of the peptide up to the second linked position.
If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos.
Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
a precursor mass for the two peptides and the linker is needed.
In the case of a loop link the precursor mass is the mass of the only peptide and the linker.
Although this function is more general, currently it is mainly used for loop-links and mono-links,
because residues in the second, unknown peptide cannot be considered for possible neutral losses.
The generated ion types and other additional settings are determined by the tool parameters
:param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
:param peptide: The peptide to fragment
:param link_pos: The position of the cross-linker on the given peptide
:param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
:param frag_alpha: True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
:param mincharge: The minimal charge of the ions
:param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
:param link_pos_2: A second position for the linker, in case it is a loop link
)doc")
        ;

    // -----------------------------------------------------------------------
    // ThresholdMower
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ThresholdMower, OpenMS::DefaultParamHandler>(m, "ThresholdMower", 
        R"doc(
Removes all peaks below an intensity threshold
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ThresholdMower &>())
        .def("__copy__", [](const OpenMS::ThresholdMower& self) { return OpenMS::ThresholdMower(self); })
        .def("__deepcopy__", [](const OpenMS::ThresholdMower& self, nb::dict) { return OpenMS::ThresholdMower(self); }, "memo"_a)
        .def("filterPeakSpectrum", [](OpenMS::ThresholdMower& self, OpenMS::MSSpectrum& spectrum) { return self.filterPeakSpectrum(spectrum); }, "spectrum"_a)
        .def("filterPeakMap", [](OpenMS::ThresholdMower& self, OpenMS::MSExperiment& exp) { return self.filterPeakMap(exp); }, "exp"_a)
        ;

    // -----------------------------------------------------------------------
    // TransitionTSVFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransitionTSVFile, OpenMS::ProgressLogger>(m, "TransitionTSVFile", 
        R"doc(
This class supports reading and writing of OpenSWATH transition lists
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("validateTargetedExperiment", [](OpenMS::TransitionTSVFile& self, const OpenMS::TargetedExperiment& targeted_exp) { return self.validateTargetedExperiment(targeted_exp); }, "targeted_exp"_a, "Validate a TargetedExperiment (check that all ids are unique)")
        .def("convertTargetedExperimentToTSV", [](OpenMS::TransitionTSVFile& self, const char* filename, OpenMS::TargetedExperiment& targeted_exp) { self.convertTargetedExperimentToTSV(filename, targeted_exp); }, "filename"_a, "targeted_exp"_a, "Write a TargetedExperiment to a TSV file")
        .def("convertTSVToTargetedExperiment", [](OpenMS::TransitionTSVFile& self, const char* filename, OpenMS::TargetedExperiment& targeted_exp) {
            self.convertTSVToTargetedExperiment(filename, OpenMS::FileTypes::TSV, targeted_exp);
        }, "filename"_a, "targeted_exp"_a, "Read a TSV file into a TargetedExperiment")
        ;

    // -----------------------------------------------------------------------
    // TransitionPQPFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransitionPQPFile, OpenMS::TransitionTSVFile>(m, "TransitionPQPFile", "This class supports reading and writing of PQP files")
        .def(nb::init<>())
        .def("validateTargetedExperiment", [](OpenMS::TransitionPQPFile& self, const OpenMS::TargetedExperiment& targeted_exp) { return self.validateTargetedExperiment(targeted_exp); }, "targeted_exp"_a)
        .def("convertPQPToTargetedExperiment", [](OpenMS::TransitionPQPFile& self, const char* filename, OpenMS::TargetedExperiment& targeted_exp, bool legacy_traml_id) { self.convertPQPToTargetedExperiment(filename, targeted_exp, legacy_traml_id); }, "filename"_a, "targeted_exp"_a, "legacy_traml_id"_a = false)
        .def("convertTargetedExperimentToPQP", [](OpenMS::TransitionPQPFile& self, const char* filename, OpenMS::TargetedExperiment& targeted_exp) { self.convertTargetedExperimentToPQP(filename, targeted_exp); }, "filename"_a, "targeted_exp"_a)
        ;

    // -----------------------------------------------------------------------
    // WindowMower
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::WindowMower, OpenMS::DefaultParamHandler>(m, "WindowMower", 
        R"doc(
Retains the highest peaks in a sliding or jumping window
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::WindowMower &>())
        .def("__copy__", [](const OpenMS::WindowMower& self) { return OpenMS::WindowMower(self); })
        .def("__deepcopy__", [](const OpenMS::WindowMower& self, nb::dict) { return OpenMS::WindowMower(self); }, "memo"_a)
        .def("filterPeakSpectrum", [](OpenMS::WindowMower& self, OpenMS::MSSpectrum& spectrum) { return self.filterPeakSpectrum(spectrum); }, "spectrum"_a)
        .def("filterPeakMap", [](OpenMS::WindowMower& self, OpenMS::MSExperiment& exp) { return self.filterPeakMap(exp); }, "exp"_a)

        .def("filterPeakSpectrumForTopNInSlidingWindow", [](OpenMS::WindowMower& self, OpenMS::MSSpectrum& spectrum) {
            self.filterPeakSpectrumForTopNInSlidingWindow(spectrum);
        }, "spectrum"_a, "Sliding window version (slower)")
        .def("filterPeakSpectrumForTopNInJumpingWindow", [](OpenMS::WindowMower& self, OpenMS::MSSpectrum& spectrum) { self.filterPeakSpectrumForTopNInJumpingWindow(spectrum); }, "spectrum"_a)
        ;

    // -----------------------------------------------------------------------
    // XFDRAlgorithm
    // -----------------------------------------------------------------------
    auto xfdralgorithm_class = nb::class_<OpenMS::XFDRAlgorithm, OpenMS::DefaultParamHandler>(m, "XFDRAlgorithm", 
        R"doc(
Calculates false discovery rate estimates on crosslink
identifications
DefaultParamHandler
)doc")
        .def(nb::init<>())
        .def("run", [](OpenMS::XFDRAlgorithm& self, OpenMS::PeptideIdentificationList& peptide_ids, OpenMS::ProteinIdentification& protein_id) { return self.run(peptide_ids, protein_id); }, "peptide_ids"_a, "protein_id"_a)
        .def("validateClassArguments", [](const OpenMS::XFDRAlgorithm& self) { return self.validateClassArguments(); })
        ;
    // ExitCodes enum nested under XFDRAlgorithm
    nb::enum_<OpenMS::XFDRAlgorithm::ExitCodes>(xfdralgorithm_class, "ExitCodes", nb::is_arithmetic())
        .value("EXECUTION_OK", OpenMS::XFDRAlgorithm::ExitCodes::EXECUTION_OK)
        .value("ILLEGAL_PARAMETERS", OpenMS::XFDRAlgorithm::ExitCodes::ILLEGAL_PARAMETERS)
        .value("UNEXPECTED_RESULT", OpenMS::XFDRAlgorithm::ExitCodes::UNEXPECTED_RESULT)
        ;

    // -----------------------------------------------------------------------
    // XMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Internal::XMLFile>(m, "XMLFile", "OpenMS class XMLFile")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::Internal::XMLFile& self) { return OpenMS::Internal::XMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::Internal::XMLFile& self, nb::dict) { return OpenMS::Internal::XMLFile(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, OpenMS::String>())
        .def("getVersion", [](const OpenMS::Internal::XMLFile& self) { return self.getVersion(); }, "Return the version of the schema")
        ;

    // -----------------------------------------------------------------------
    // ConsensusXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusXMLFile, OpenMS::Internal::XMLFile>(m, "ConsensusXMLFile", 
        R"doc(
File adapter for consensusXML files
Provides methods to load and store consensus maps in consensusXML format.
ConsensusXML files store consensus features linking corresponding features across
multiple LC-MS runs, typically used for quantification workflows.
Usage:
.. code-block:: python
cm = ConsensusMap()
ConsensusXMLFile().load("test.consensusXML", cm)
for cf in cm:
print(cf.getRT(), cf.getMZ(), cf.getIntensity())
)doc")
        .def(nb::init<>())
        .def("getOptions", [](OpenMS::ConsensusXMLFile& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Mutable access to the options for loading/storing")

        .def("load", [](OpenMS::ConsensusXMLFile& self, const OpenMS::String& filename, OpenMS::ConsensusMap& map) {
            nb::gil_scoped_release release;
            self.load(filename, map);
        }, "filename"_a, "map"_a, "Load a consensusXML file")

        .def("store", [](OpenMS::ConsensusXMLFile& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& map) {
            nb::gil_scoped_release release;
            self.store(filename, map);
        }, "filename"_a, "map"_a, "Store a ConsensusMap to consensusXML")
        ;

    // -----------------------------------------------------------------------
    // FeatureXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureXMLFile, OpenMS::Internal::XMLFile>(m, "FeatureXMLFile", 
        R"doc(
File adapter for featureXML files
Provides methods to load and store feature maps in featureXML format.
FeatureXML files store LC-MS features with their convex hulls, intensities, and metadata.
Usage:
.. code-block:: python
fm = FeatureMap()
FeatureXMLFile().load("test.featureXML", fm)
for feature in fm:
print(feature.getRT(), feature.getMZ(), feature.getIntensity())
)doc")
        .def(nb::init<>())
        .def("loadSize", [](OpenMS::FeatureXMLFile& self, const OpenMS::String& filename) { return self.loadSize(filename); }, "filename"_a, "Counts the number of features in the file without loading the full data")
        .def("getOptions", [](OpenMS::FeatureXMLFile& self) -> OpenMS::FeatureFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Access to the options for loading/storing")
        .def("setOptions", [](OpenMS::FeatureXMLFile& self, const OpenMS::FeatureFileOptions& p0) { return self.setOptions(p0); }, "Setter for options for loading/storing")

        .def("load", [](OpenMS::FeatureXMLFile& self, const OpenMS::String& filename, OpenMS::FeatureMap& map) {
            nb::gil_scoped_release release;
            self.load(filename, map);
        }, "filename"_a, "map"_a, "Load a featureXML file")

        .def("store", [](OpenMS::FeatureXMLFile& self, const OpenMS::String& filename, const OpenMS::FeatureMap& map) {
            nb::gil_scoped_release release;
            self.store(filename, map);
        }, "filename"_a, "map"_a, "Store a FeatureMap to featureXML")
        ;

    // -----------------------------------------------------------------------
    // IdXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IdXMLFile, OpenMS::Internal::XMLFile>(m, "IdXMLFile", 
        R"doc(
File adapter for idXML files
Provides methods to load and store identification data in idXML format.
idXML files store protein and peptide identifications from database search engines.
Usage:
.. code-block:: python
protein_ids = []
peptide_ids = []
IdXMLFile().load("test.idXML", protein_ids, peptide_ids)
)doc")
        .def(nb::init<>())

        .def("_load_internal", [](OpenMS::IdXMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            {
                nb::gil_scoped_release release;
                self.load(filename, proteins, peptides);
            }
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an idXML file, returns tuple (proteins, peptides)")

        .def("_store_internal", [](OpenMS::IdXMLFile& self, const OpenMS::String& filename,
                         const std::vector<OpenMS::ProteinIdentification>& proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            nb::gil_scoped_release release;
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to an idXML file")
        ;

    // -----------------------------------------------------------------------
    // MascotXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MascotXMLFile, OpenMS::Internal::XMLFile>(m, "MascotXMLFile", 
        R"doc(
Used to load Mascot XML files
XMLFile
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MascotXMLFile& self) { return OpenMS::MascotXMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::MascotXMLFile& self, nb::dict) { return OpenMS::MascotXMLFile(self); }, "memo"_a)
        .def_static("initializeLookup", [](OpenMS::SpectrumMetaDataLookup& lookup, const OpenMS::PeakMap& experiment, const OpenMS::String& scan_regex) { OpenMS::MascotXMLFile::initializeLookup(lookup, experiment, scan_regex); }, "lookup"_a, "experiment"_a, "scan_regex"_a = "", "Initialize spectrum lookup")
        .def("load", [](OpenMS::MascotXMLFile& self, const OpenMS::String& filename, OpenMS::ProteinIdentification& protein_identification, OpenMS::PeptideIdentificationList& id_data, const OpenMS::SpectrumMetaDataLookup& lookup) { self.load(filename, protein_identification, id_data, lookup); }, "filename"_a, "protein_identification"_a, "id_data"_a, "lookup"_a, "Load Mascot XML file")
        ;

    // -----------------------------------------------------------------------
    // MzDataFile
    // -----------------------------------------------------------------------
    auto mzdatafile_class = nb::class_<OpenMS::MzDataFile, OpenMS::Internal::XMLFile>(m, "MzDataFile", 
        R"doc(
File adapter for MzData files
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzDataFile& self) { return OpenMS::MzDataFile(self); })
        .def("__deepcopy__", [](const OpenMS::MzDataFile& self, nb::dict) { return OpenMS::MzDataFile(self); }, "memo"_a)
        .def("getOptions", [](OpenMS::MzDataFile& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Returns the options for loading/storing")
        .def("setOptions", [](OpenMS::MzDataFile& self, const OpenMS::PeakFileOptions& p0) { return self.setOptions(p0); }, "Sets options for loading/storing")
        .def("load", [](OpenMS::MzDataFile& self, const OpenMS::String& filename) { OpenMS::MSExperiment map; { nb::gil_scoped_release release; self.load(filename, map); } return map; }, "filename"_a)
        .def("load", [](OpenMS::MzDataFile& self, const OpenMS::String& filename, OpenMS::MSExperiment& map) { nb::gil_scoped_release release; self.load(filename, map); }, "filename"_a, "map"_a, "Loads a map from a MzData file into the given MSExperiment")
        .def("store", [](const OpenMS::MzDataFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& map) { nb::gil_scoped_release release; return self.store(filename, map); }, "filename"_a, "map"_a,
            R"doc(
Loads a map from a MzData file
:param filename: Directory of the file with the file name
:param map: It has to be a MSExperiment or have the same interface
:raises:
Exception: FileNotFound is thrown if the file could not be opened
:raises:
Exception: ParseError is thrown if an error occurs during parsing
)doc")
        .def("isSemanticallyValid", [](OpenMS::MzDataFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::String> errors;
            std::vector<OpenMS::String> warnings;
            auto result = self.isSemanticallyValid(filename, errors, warnings);
            return nb::make_tuple(result, errors, warnings);
        }, "filename"_a)
        ;
    def_ProgressLogger<OpenMS::MzDataFile>(mzdatafile_class);

    // -----------------------------------------------------------------------
    // MzIdentMLFile
    // -----------------------------------------------------------------------
    auto mzidentmlfile_class = nb::class_<OpenMS::MzIdentMLFile, OpenMS::Internal::XMLFile>(m, "MzIdentMLFile", 
        R"doc(
File adapter for MzIdentML files
ProgressLogger
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzIdentMLFile& self) { return OpenMS::MzIdentMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::MzIdentMLFile& self, nb::dict) { return OpenMS::MzIdentMLFile(self); }, "memo"_a)
        .def("isSemanticallyValid", [](OpenMS::MzIdentMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::String> errors;
            std::vector<OpenMS::String> warnings;
            auto result = self.isSemanticallyValid(filename, errors, warnings);
            return nb::make_tuple(result, errors, warnings);
        }, "filename"_a, 
            R"doc(
Stores the identifications in a MzIdentML file
:raises:
Exception: UnableToCreateFile is thrown if the file could not be created
)doc")

        .def("_load_internal", [](OpenMS::MzIdentMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an mzIdentML file, returns tuple (proteins, peptides)")

        .def("_store_internal", [](OpenMS::MzIdentMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to an mzIdentML file")
        ;
    def_ProgressLogger<OpenMS::MzIdentMLFile>(mzidentmlfile_class);

    // -----------------------------------------------------------------------
    // MzMLFile
    // -----------------------------------------------------------------------
    auto mzmlfile_class = nb::class_<OpenMS::MzMLFile, OpenMS::Internal::XMLFile>(m, "MzMLFile", 
        R"doc(
ProgressLogger

File adapter for MzML files
Provides methods to load and store MzML files.
PeakFileOptions allow to load a reduced subset of the data into an MSExperiment.
See help(MSExperiment) how data is stored after loading.
See help(PeakFileOptions) for available options.
Usage:
.. code-block:: python
exp = MSExperiment()
MzMLFile().load("test.mzML", exp)
spec = []
for s in exp.getSpectra():
if s.getMSLevel() != 1:
spec.append(s)
exp.setSpectra(spec)
MzMLFile().store("filtered.mzML", exp)
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzMLFile& self) { return OpenMS::MzMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::MzMLFile& self, nb::dict) { return OpenMS::MzMLFile(self); }, "memo"_a)
        .def("getOptions", [](OpenMS::MzMLFile& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Returns the options for loading/storing")
        .def("setOptions", [](OpenMS::MzMLFile& self, const OpenMS::PeakFileOptions& p0) { return self.setOptions(p0); }, "Set PeakFileOptions to perform filtering during loading. E.g., to load only MS1 spectra or meta data only")

        .def("load", [](OpenMS::MzMLFile& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) {
            nb::gil_scoped_release release;
            self.load(filename, exp);
        }, "filename"_a, "exp"_a, "Load an mzML file into an MSExperiment")

        .def("store", [](OpenMS::MzMLFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp) {
            nb::gil_scoped_release release;
            self.store(filename, exp);
        }, "filename"_a, "exp"_a, "Store an MSExperiment to an mzML file")

        .def("storeBuffer", [](OpenMS::MzMLFile& self, nb::object output_str, const OpenMS::MSExperiment& exp) {
            std::string buf;
            {
                nb::gil_scoped_release release;
                self.storeBuffer(buf, exp);
            }
            output_str.attr("_value") = nb::cast(buf);
        }, "output"_a, "exp"_a, "Store an MSExperiment to an in-memory mzML string buffer")

        .def("loadBuffer", [](OpenMS::MzMLFile& self, const std::string& buffer, OpenMS::MSExperiment& exp) {
            nb::gil_scoped_release release;
            self.loadBuffer(buffer, exp);
        }, "buffer"_a, "exp"_a, "Load an MSExperiment from an in-memory mzML string buffer")

        .def("transform", [](OpenMS::MzMLFile& self, const OpenMS::String& filename, nb::object consumer,
                             bool skip_full_count, bool skip_first_pass) {
            NanobindMSDataConsumer wrapper(consumer);
            self.transform(filename, &wrapper, skip_full_count, skip_first_pass);
        }, "filename"_a, "consumer"_a, "skip_full_count"_a = false, "skip_first_pass"_a = false,
        "Transform an mzML file using a consumer object (streaming processing)")
        .def("isSemanticallyValid", [](OpenMS::MzMLFile& self, const OpenMS::String& filename) { OpenMS::StringList errors; OpenMS::StringList warnings; bool result = self.isSemanticallyValid(filename, errors, warnings); return nb::make_tuple(result, errors, warnings); }, "filename"_a, "Check semantic validity and return (is_valid, errors, warnings)")
        ;
    def_ProgressLogger<OpenMS::MzMLFile>(mzmlfile_class);

    // -----------------------------------------------------------------------
    // MzXMLFile
    // -----------------------------------------------------------------------
    auto mzxmlfile_class = nb::class_<OpenMS::MzXMLFile, OpenMS::Internal::XMLFile>(m, "MzXMLFile", 
        R"doc(
ProgressLogger

File adapter for MzXML files
Provides methods to load and store MzXML files.
MzXML is an older format; for new projects consider using MzML instead.
Usage:
.. code-block:: python
exp = MSExperiment()
MzXMLFile().load("test.mzXML", exp)
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzXMLFile& self) { return OpenMS::MzXMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::MzXMLFile& self, nb::dict) { return OpenMS::MzXMLFile(self); }, "memo"_a)
        .def("getOptions", [](OpenMS::MzXMLFile& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Returns the options for loading/storing")
        .def("setOptions", [](OpenMS::MzXMLFile& self, const OpenMS::PeakFileOptions& p0) { return self.setOptions(p0); }, "Sets options for loading/storing")

        .def("load", [](OpenMS::MzXMLFile& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) {
            nb::gil_scoped_release release;
            self.load(filename, exp);
        }, "filename"_a, "exp"_a, "Load an mzXML file into an MSExperiment")

        .def("store", [](OpenMS::MzXMLFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp) {
            nb::gil_scoped_release release;
            self.store(filename, exp);
        }, "filename"_a, "exp"_a, "Store an MSExperiment to an mzXML file")

        .def("transform", [](OpenMS::MzXMLFile& self, const OpenMS::String& filename, nb::object consumer,
                             bool skip_full_count) {
            NanobindMSDataConsumer wrapper(consumer);
            self.transform(filename, &wrapper, skip_full_count);
        }, "filename"_a, "consumer"_a, "skip_full_count"_a = false,
        "Transform an mzXML file using a consumer object (streaming processing)")
        ;
    def_ProgressLogger<OpenMS::MzXMLFile>(mzxmlfile_class);

    // -----------------------------------------------------------------------
    // OMSSAXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OMSSAXMLFile, OpenMS::Internal::XMLFile>(m, "OMSSAXMLFile", 
        R"doc(
XMLFile

Used to load OMSSAXML files
This class is used to load documents that implement
the schema of OMSSAXML files
)doc")
        .def(nb::init<>())
        .def("load", [](OpenMS::OMSSAXMLFile& self, const OpenMS::String& filename, bool load_proteins, bool load_empty_hits) { OpenMS::ProteinIdentification protein_identification; OpenMS::PeptideIdentificationList id_data; self.load(filename, protein_identification, id_data, load_proteins, load_empty_hits); return std::make_tuple(protein_identification, id_data); }, "filename"_a, "load_proteins"_a, "load_empty_hits"_a)
        .def("setModificationDefinitionsSet", [](OpenMS::OMSSAXMLFile& self, const OpenMS::ModificationDefinitionsSet& rhs) { return self.setModificationDefinitionsSet(rhs); }, "rhs"_a, "Sets the valid modifications")
        ;

    // -----------------------------------------------------------------------
    // ParamXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ParamXMLFile, OpenMS::Internal::XMLFile>(m, "ParamXMLFile", 
        R"doc(
The file pendant of the Param class used to load and store the param
datastructure as paramXML
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::ParamXMLFile& self) { return OpenMS::ParamXMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::ParamXMLFile& self, nb::dict) { return OpenMS::ParamXMLFile(self); }, "memo"_a)
        .def("store", [](const OpenMS::ParamXMLFile& self, const OpenMS::String& filename, const OpenMS::Param& param) { return self.store(filename, param); }, "filename"_a, "param"_a, 
            R"doc(
Read XML file
:param filename: The file from where to read the Param object
:param param: The param object where the read data should be stored
:raises:
Exception: FileNotFound is thrown if the file could not be found
:raises:
Exception: ParseError is thrown if an error occurs during parsing
)doc")
        .def("load", [](OpenMS::ParamXMLFile& self, const OpenMS::String& filename) { OpenMS::Param param; self.load(filename, param); return param; }, "filename"_a)
        .def("load", [](OpenMS::ParamXMLFile& self, const OpenMS::String& filename, OpenMS::Param& param) { self.load(filename, param); }, "filename"_a, "param"_a)
        ;

    // -----------------------------------------------------------------------
    // PepXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PepXMLFile, OpenMS::Internal::XMLFile>(m, "PepXMLFile", "Used to load and store PepXML files")
        .def(nb::init<>())
        .def("keepNativeSpectrumName", [](OpenMS::PepXMLFile& self, bool keep) { return self.keepNativeSpectrumName(keep); }, "keep"_a)
        .def("setParseUnknownScores", [](OpenMS::PepXMLFile& self, bool parse_unknown_scores) { return self.setParseUnknownScores(parse_unknown_scores); }, "parse_unknown_scores"_a)

        .def("_load_internal", [](OpenMS::PepXMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load a pepXML file, returns tuple (proteins, peptides)")

        .def("_store_internal", [](OpenMS::PepXMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to a pepXML file")
        ;

    // -----------------------------------------------------------------------
    // PepXMLFileMascot
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PepXMLFileMascot, OpenMS::Internal::XMLFile>(m, "PepXMLFileMascot", 
        R"doc(
Used to load Mascot PepXML files
A schema for this format can be found at http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd
)doc")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // ProtXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProtXMLFile, OpenMS::Internal::XMLFile>(m, "ProtXMLFile", 
        R"doc(
Used to load (storing not supported, yet) ProtXML files
This class is used to load (storing not supported, yet) documents that implement
the schema of ProtXML files
A documented schema for this format comes with the TPP and can also be
found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS
OpenMS can only read parts of the protein_summary subtree to extract
protein-peptide associations. All other parts are silently ignored
For protein groups, only the "group leader" (which is annotated with a
probability and coverage) receives these attributes. All indistinguishable
proteins of the same group only have an accession and score of -1
)doc")
        .def(nb::init<>())
        .def("load", [](OpenMS::ProtXMLFile& self, const OpenMS::String& filename) { OpenMS::ProteinIdentification protein_ids; OpenMS::PeptideIdentification peptide_ids; self.load(filename, protein_ids, peptide_ids); return std::make_tuple(protein_ids, peptide_ids); }, "filename"_a)
        .def("store", [](OpenMS::ProtXMLFile& self, const OpenMS::String& filename, const OpenMS::ProteinIdentification& protein_ids, const OpenMS::PeptideIdentification& peptide_ids, const OpenMS::String& document_id) { return self.store(filename, protein_ids, peptide_ids, document_id); }, "filename"_a, "protein_ids"_a, "peptide_ids"_a, "document_id"_a = "", 
            R"doc(
Loads the identifications of an ProtXML file without identifier
The information is read in and the information is stored in the
corresponding variables
:raises:
Exception: FileNotFound is thrown if the file could not be found
:raises:
Exception: ParseError is thrown if an error occurs during parsing
Not implemented
)doc")
        ;

    // -----------------------------------------------------------------------
    // QualityParameter (QcMLFile::QualityParameter)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::QcMLFile::QualityParameter>(m, "QualityParameter",
        "Representation of a quality parameter in QcMLFile")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::QcMLFile::QualityParameter &>())
        .def("__copy__", [](const OpenMS::QcMLFile::QualityParameter& self) { return OpenMS::QcMLFile::QualityParameter(self); })
        .def("__deepcopy__", [](const OpenMS::QcMLFile::QualityParameter& self, nb::dict) { return OpenMS::QcMLFile::QualityParameter(self); }, "memo"_a)
        .def_rw("name", &OpenMS::QcMLFile::QualityParameter::name)
        .def_rw("id", &OpenMS::QcMLFile::QualityParameter::id)
        .def_rw("value", &OpenMS::QcMLFile::QualityParameter::value)
        .def_rw("cvRef", &OpenMS::QcMLFile::QualityParameter::cvRef)
        .def_rw("cvAcc", &OpenMS::QcMLFile::QualityParameter::cvAcc)
        .def_rw("unitRef", &OpenMS::QcMLFile::QualityParameter::unitRef)
        .def_rw("unitAcc", &OpenMS::QcMLFile::QualityParameter::unitAcc)
        .def_rw("flag", &OpenMS::QcMLFile::QualityParameter::flag)
        .def("toXMLString", [](const OpenMS::QcMLFile::QualityParameter& self, OpenMS::UInt indentation_level) { return self.toXMLString(indentation_level); }, "indentation_level"_a)
        .def(nb::self == nb::self)
        .def(nb::self < nb::self)
        ;

    // -----------------------------------------------------------------------
    // Attachment (QcMLFile::Attachment)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::QcMLFile::Attachment>(m, "Attachment",
        "Attachment for QcMLFile quality parameters")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::QcMLFile::Attachment &>())
        .def("__copy__", [](const OpenMS::QcMLFile::Attachment& self) { return OpenMS::QcMLFile::Attachment(self); })
        .def("__deepcopy__", [](const OpenMS::QcMLFile::Attachment& self, nb::dict) { return OpenMS::QcMLFile::Attachment(self); }, "memo"_a)
        .def_rw("name", &OpenMS::QcMLFile::Attachment::name)
        .def_rw("id", &OpenMS::QcMLFile::Attachment::id)
        .def_rw("value", &OpenMS::QcMLFile::Attachment::value)
        .def_rw("cvRef", &OpenMS::QcMLFile::Attachment::cvRef)
        .def_rw("cvAcc", &OpenMS::QcMLFile::Attachment::cvAcc)
        .def_rw("unitRef", &OpenMS::QcMLFile::Attachment::unitRef)
        .def_rw("unitAcc", &OpenMS::QcMLFile::Attachment::unitAcc)
        .def_rw("binary", &OpenMS::QcMLFile::Attachment::binary)
        .def_rw("qualityRef", &OpenMS::QcMLFile::Attachment::qualityRef)
        .def_rw("colTypes", &OpenMS::QcMLFile::Attachment::colTypes)
        .def_rw("tableRows", &OpenMS::QcMLFile::Attachment::tableRows)
        .def("toXMLString", [](const OpenMS::QcMLFile::Attachment& self, OpenMS::UInt indentation_level) { return self.toXMLString(indentation_level); }, "indentation_level"_a)
        .def("toCSVString", [](const OpenMS::QcMLFile::Attachment& self, const OpenMS::String& separator) { return self.toCSVString(separator); }, "separator"_a)
        .def(nb::self == nb::self)
        .def(nb::self < nb::self)
        ;

    // -----------------------------------------------------------------------
    // QcMLFile
    // -----------------------------------------------------------------------
    auto qcmlfile_class = nb::class_<OpenMS::QcMLFile, OpenMS::Internal::XMLFile>(m, "QcMLFile",
        R"doc(
XMLHandler
XMLFile
ProgressLogger

File adapter for QcML files used to load and store QcML files
This Class is supposed to internally collect the data for the qcML File
)doc")
        .def(nb::init<>())
        .def("map2csv", [](const OpenMS::QcMLFile& self, const std::map<OpenMS::String, std::map<OpenMS::String, OpenMS::String>>& cvs_table, const OpenMS::String& separator) { return self.map2csv(cvs_table, separator); }, "cvs_table"_a, "separator"_a, "Converts a map of maps to a CSV string")
        .def("exportIDstats", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename) { return self.exportIDstats(filename); }, "filename"_a)
        .def("registerRun", [](OpenMS::QcMLFile& self, const OpenMS::String& id, const OpenMS::String& name) { return self.registerRun(id, name); }, "id"_a, "name"_a, "Registers a run in the qcml file with the respective mappings")
        .def("registerSet", [](OpenMS::QcMLFile& self, const OpenMS::String& id, const OpenMS::String& name, const std::set<OpenMS::String>& names) { return self.registerSet(id, name, names); }, "id"_a, "name"_a, "names"_a, "Registers a set in the qcml file with the respective mappings")
        .def("addRunQualityParameter", [](OpenMS::QcMLFile& self, const OpenMS::String& r, const OpenMS::QcMLFile::QualityParameter& qp) { return self.addRunQualityParameter(r, qp); }, "r"_a, "qp"_a, "Adds a QualityParameter to run by the name r")
        .def("addRunAttachment", [](OpenMS::QcMLFile& self, const OpenMS::String& r, const OpenMS::QcMLFile::Attachment& at) { return self.addRunAttachment(r, at); }, "r"_a, "at"_a, "Adds a attachment to run by the name r")
        .def("addSetQualityParameter", [](OpenMS::QcMLFile& self, const OpenMS::String& r, const OpenMS::QcMLFile::QualityParameter& qp) { return self.addSetQualityParameter(r, qp); }, "r"_a, "qp"_a, "Adds a QualityParameter to set by the name r")
        .def("addSetAttachment", [](OpenMS::QcMLFile& self, const OpenMS::String& r, const OpenMS::QcMLFile::Attachment& at) { return self.addSetAttachment(r, at); }, "r"_a, "at"_a, "Adds a attachment to set by the name r")
        .def("removeAttachment", [](OpenMS::QcMLFile& self, const OpenMS::String& r, std::vector<OpenMS::String> ids, const OpenMS::String& at) { return self.removeAttachment(r, ids, at); }, "r"_a, "ids"_a, "at"_a = "", "Removes attachments referencing an id given in ids, from run/set r. All attachments if no attachment name is given with at")
        .def("removeAttachment", [](OpenMS::QcMLFile& self, const OpenMS::String& r, const OpenMS::String& at) { return self.removeAttachment(r, at); }, "r"_a, "at"_a, "Removes attachment with cv accession at from run/set r")
        .def("removeAllAttachments", [](OpenMS::QcMLFile& self, const OpenMS::String& at) { return self.removeAllAttachments(at); }, "at"_a, "Removes attachment with cv accession at from  all runs/sets")
        .def("removeQualityParameter", [](OpenMS::QcMLFile& self, const OpenMS::String& r, std::vector<OpenMS::String> ids) { return self.removeQualityParameter(r, ids); }, "r"_a, "ids"_a, "Removes QualityParameter going by one of the ID attributes given in ids")
        .def("merge", [](OpenMS::QcMLFile& self, const OpenMS::QcMLFile& addendum, const OpenMS::String& setname) { return self.merge(addendum, setname); }, "addendum"_a, "setname"_a = "", "Merges the given QCFile into this one")
        .def("collectSetParameter", [](OpenMS::QcMLFile& self, const OpenMS::String& setname, const OpenMS::String& qp) { std::vector<OpenMS::String> ret; self.collectSetParameter(setname, qp, ret); return ret; }, "setname"_a, "qp"_a, "Collects the values of given QPs (as CVid) of the given set")
        .def("exportAttachment", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename, const OpenMS::String& qpname) { return self.exportAttachment(filename, qpname); }, "filename"_a, "qpname"_a, "Returns a String of a tab separated rows if found empty string else from run/set by the name filename of the qualityparameter by the name qpname")
        .def("exportQP", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename, const OpenMS::String& qpname) { return self.exportQP(filename, qpname); }, "filename"_a, "qpname"_a, "Returns a String value in quotation of a QualityParameter by the name qpname in run/set by the name filename")
        .def("exportQPs", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename, const std::vector<OpenMS::String>& qpnames) { return self.exportQPs(filename, qpnames); }, "filename"_a, "qpnames"_a, "Returns a String of a tab separated QualityParameter by the name qpname in run/set by the name filename")
        .def("getRunIDs", [](const OpenMS::QcMLFile& self) { std::vector<OpenMS::String> ids; self.getRunIDs(ids); return ids; }, "Gives the ids of the registered runs in the vector ids")
        .def("getRunNames", [](const OpenMS::QcMLFile& self) { std::vector<OpenMS::String> ids; self.getRunNames(ids); return ids; }, "Gives the names of the registered runs in the vector ids")
        .def("existsRun", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename, bool checkname) { return self.existsRun(filename, checkname); }, "filename"_a, "checkname"_a = false, "Returns true if the given run id is present in this file, if checkname is true it also checks the names")
        .def("existsSet", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename, bool checkname) { return self.existsSet(filename, checkname); }, "filename"_a, "checkname"_a = false, "Returns true if the given set id is present in this file, if checkname is true it also checks the names")
        .def("existsRunQualityParameter", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename, const OpenMS::String& qpname) {
            std::vector<OpenMS::String> ids;
            self.existsRunQualityParameter(filename, qpname, ids);
            return ids;
        }, "filename"_a, "qpname"_a, "Returns the ids of the parameter name given if found in given run empty else")
        .def("existsSetQualityParameter", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename, const OpenMS::String& qpname) {
            std::vector<OpenMS::String> ids;
            self.existsSetQualityParameter(filename, qpname, ids);
            return ids;
        }, "filename"_a, "qpname"_a, "Returns the ids of the parameter name given if found in given set, empty else")
        .def("store", [](const OpenMS::QcMLFile& self, const OpenMS::String& filename) { return self.store(filename); }, "filename"_a, "Store the qcML file")
        .def("load", [](OpenMS::QcMLFile& self, const OpenMS::String& filename) { return self.load(filename); }, "filename"_a, "Load a QCFile")
        .def("reset", [](OpenMS::QcMLFile& self) { return self.reset(); })
        .def("error", [](OpenMS::QcMLFile& self, const xercesc::SAXParseException& exception) { return self.error(exception); }, "exception"_a)
        .def("warning", [](OpenMS::QcMLFile& self, const xercesc::SAXParseException& exception) { return self.warning(exception); }, "exception"_a)
        ;
    def_ProgressLogger<OpenMS::QcMLFile>(qcmlfile_class);

    // -----------------------------------------------------------------------
    // TraMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TraMLFile, OpenMS::Internal::XMLFile>(m, "TraMLFile", 
        R"doc(
File adapter for TraML files
Provides methods to load and store targeted experiment data in TraML format.
TraML (Transition Markup Language) files store SRM/MRM transitions for targeted
proteomics and metabolomics experiments.
Usage:
.. code-block:: python
targeted_exp = TargetedExperiment()
TraMLFile().load("transitions.TraML", targeted_exp)
for transition in targeted_exp.getTransitions():
print(transition.getPrecursorMZ(), transition.getProductMZ())
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::TraMLFile& self) { return OpenMS::TraMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::TraMLFile& self, nb::dict) { return OpenMS::TraMLFile(self); }, "memo"_a)

        .def("load", [](OpenMS::TraMLFile& self, const OpenMS::String& filename, OpenMS::TargetedExperiment& exp) {
            nb::gil_scoped_release release;
            self.load(filename, exp);
        }, "filename"_a, "exp"_a, "Load a TraML file")

        .def("store", [](OpenMS::TraMLFile& self, const OpenMS::String& filename, const OpenMS::TargetedExperiment& exp) {
            nb::gil_scoped_release release;
            self.store(filename, exp);
        }, "filename"_a, "exp"_a, "Store to a TraML file")
        .def("isSemanticallyValid", [](OpenMS::TraMLFile& self, const OpenMS::String& filename) { OpenMS::StringList errors; OpenMS::StringList warnings; bool result = self.isSemanticallyValid(filename, errors, warnings); return nb::make_tuple(result, errors, warnings); }, "filename"_a, "Check semantic validity and return (is_valid, errors, warnings)")
        ;

    // -----------------------------------------------------------------------
    // TransformationXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationXMLFile, OpenMS::Internal::XMLFile>(m, "TransformationXMLFile", "Used to load and store TransformationXML files")
        .def(nb::init<>())
        .def("load", [](OpenMS::TransformationXMLFile& self, const OpenMS::String& filename, bool fit_model) { OpenMS::TransformationDescription transformation; self.load(filename, transformation, fit_model); return transformation; }, "filename"_a, "fit_model"_a)
        .def("load", [](OpenMS::TransformationXMLFile& self, const OpenMS::String& filename, OpenMS::TransformationDescription& td, bool fit_model) { self.load(filename, td, fit_model); }, "filename"_a, "td"_a, "fit_model"_a)
        .def("store", [](OpenMS::TransformationXMLFile& self, const OpenMS::String& filename, const OpenMS::TransformationDescription& transformation) { return self.store(filename, transformation); }, "filename"_a, "transformation"_a)
        ;

    // -----------------------------------------------------------------------
    // UnimodXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::UnimodXMLFile, OpenMS::Internal::XMLFile>(m, "UnimodXMLFile", 
        R"doc(
Used to load XML files from unimod.org files
XMLFile
)doc")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // XQuestResultXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::XQuestResultXMLFile, OpenMS::Internal::XMLFile>(m, "XQuestResultXMLFile", 
        R"doc(
Used to load and store xQuest result files
XMLFile
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::XQuestResultXMLFile& self) { return OpenMS::XQuestResultXMLFile(self); })
        .def("__deepcopy__", [](const OpenMS::XQuestResultXMLFile& self, nb::dict) { return OpenMS::XQuestResultXMLFile(self); }, "memo"_a)
        .def("getNumberOfHits", [](const OpenMS::XQuestResultXMLFile& self) { return self.getNumberOfHits(); }, "Returns the total number of hits in the file")
        .def("getMinScore", [](const OpenMS::XQuestResultXMLFile& self) { return self.getMinScore(); }, "Returns minimum score among the hits in the file")
        .def("getMaxScore", [](const OpenMS::XQuestResultXMLFile& self) { return self.getMaxScore(); }, "Returns maximum score among the hits in the file")

        .def("_load_internal", [](OpenMS::XQuestResultXMLFile& self, const OpenMS::String& filename) {
            OpenMS::PeptideIdentificationList peptides;
            std::vector<OpenMS::ProteinIdentification> proteins;
            self.load(filename, peptides, proteins);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an xQuest result XML file, returns tuple (proteins, peptides)")

        .def("_store_internal", [](const OpenMS::XQuestResultXMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to an xQuest result XML file")
        .def("load", [](OpenMS::XQuestResultXMLFile& self, const OpenMS::String& filename, OpenMS::PeptideIdentificationList& pep_ids, std::vector<OpenMS::ProteinIdentification>& prot_ids) { self.load(filename, pep_ids, prot_ids); }, "filename"_a, "pep_ids"_a, "prot_ids"_a)
        .def("store", [](const OpenMS::XQuestResultXMLFile& self, const OpenMS::String& filename, const std::vector<OpenMS::ProteinIdentification>& poid, const OpenMS::PeptideIdentificationList& peid) { self.store(filename, poid, peid); }, "filename"_a, "poid"_a, "peid"_a)
        .def_static("writeXQuestXMLSpec", [](
                const OpenMS::String& out_file,
                const OpenMS::String& base_name,
                const OpenMS::OPXLDataStructs::PreprocessedPairSpectra& preprocessed_pair_spectra,
                const std::vector<std::pair<size_t, size_t>>& spectrum_pairs,
                const std::vector<std::vector<OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch>>& all_top_csms,
                const OpenMS::PeakMap& spectra,
                bool test_mode) {
            OpenMS::XQuestResultXMLFile::writeXQuestXMLSpec(out_file, base_name,
                preprocessed_pair_spectra, spectrum_pairs, all_top_csms, spectra, test_mode);
        }, "out_file"_a, "base_name"_a, "preprocessed_pair_spectra"_a, "spectrum_pairs"_a,
           "all_top_csms"_a, "spectra"_a, "test_mode"_a = false,
           "Write xQuest XML spec file (labeled, with pair spectra)")
        .def_static("writeXQuestXMLSpec", [](
                const OpenMS::String& out_file,
                const OpenMS::String& base_name,
                const std::vector<std::vector<OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch>>& all_top_csms,
                const OpenMS::PeakMap& spectra,
                bool test_mode) {
            OpenMS::XQuestResultXMLFile::writeXQuestXMLSpec(out_file, base_name,
                all_top_csms, spectra, test_mode);
        }, "out_file"_a, "base_name"_a, "all_top_csms"_a, "spectra"_a, "test_mode"_a = false,
           "Write xQuest XML spec file (label-free)")
        ;

    // -----------------------------------------------------------------------
    // XTandemInfile
    // -----------------------------------------------------------------------
    auto xtandeminfile_class = nb::class_<OpenMS::XTandemInfile, OpenMS::Internal::XMLFile>(m, "XTandemInfile", "XTandem input file")
        .def(nb::init<>())
        .def("setFragmentMassTolerance", [](OpenMS::XTandemInfile& self, double tolerance) { return self.setFragmentMassTolerance(tolerance); }, "tolerance"_a)
        .def("getFragmentMassTolerance", [](const OpenMS::XTandemInfile& self) { return self.getFragmentMassTolerance(); })
        .def("setPrecursorMassTolerancePlus", [](OpenMS::XTandemInfile& self, double tol) { return self.setPrecursorMassTolerancePlus(tol); }, "tol"_a)
        .def("getPrecursorMassTolerancePlus", [](const OpenMS::XTandemInfile& self) { return self.getPrecursorMassTolerancePlus(); })
        .def("setPrecursorMassToleranceMinus", [](OpenMS::XTandemInfile& self, double tol) { return self.setPrecursorMassToleranceMinus(tol); }, "tol"_a)
        .def("getPrecursorMassToleranceMinus", [](const OpenMS::XTandemInfile& self) { return self.getPrecursorMassToleranceMinus(); })
        .def("setPrecursorErrorType", [](OpenMS::XTandemInfile& self, OpenMS::XTandemInfile::MassType mono_isotopic) { return self.setPrecursorErrorType(mono_isotopic); }, "mono_isotopic"_a)
        .def("getPrecursorErrorType", [](const OpenMS::XTandemInfile& self) { return self.getPrecursorErrorType(); })
        .def("setFragmentMassErrorUnit", [](OpenMS::XTandemInfile& self, OpenMS::XTandemInfile::ErrorUnit unit) { return self.setFragmentMassErrorUnit(unit); }, "unit"_a)
        .def("getFragmentMassErrorUnit", [](const OpenMS::XTandemInfile& self) { return self.getFragmentMassErrorUnit(); })
        .def("setPrecursorMassErrorUnit", [](OpenMS::XTandemInfile& self, OpenMS::XTandemInfile::ErrorUnit unit) { return self.setPrecursorMassErrorUnit(unit); }, "unit"_a)
        .def("getPrecursorMassErrorUnit", [](const OpenMS::XTandemInfile& self) { return self.getPrecursorMassErrorUnit(); })
        .def("setNumberOfThreads", [](OpenMS::XTandemInfile& self, unsigned int threads) { return self.setNumberOfThreads(threads); }, "threads"_a)
        .def("getNumberOfThreads", [](const OpenMS::XTandemInfile& self) { return self.getNumberOfThreads(); })
        .def("setModifications", [](OpenMS::XTandemInfile& self, const OpenMS::ModificationDefinitionsSet& mods) { return self.setModifications(mods); }, "mods"_a)
        .def("getModifications", [](const OpenMS::XTandemInfile& self) -> const OpenMS::ModificationDefinitionsSet & { return self.getModifications(); }, nb::rv_policy::reference_internal)
        .def("setOutputFilename", [](OpenMS::XTandemInfile& self, const OpenMS::String& output) { return self.setOutputFilename(output); }, "output"_a)
        .def("getOutputFilename", [](const OpenMS::XTandemInfile& self) { return self.getOutputFilename(); })
        .def("setInputFilename", [](OpenMS::XTandemInfile& self, const OpenMS::String& input_file) { return self.setInputFilename(input_file); }, "input_file"_a)
        .def("getInputFilename", [](const OpenMS::XTandemInfile& self) { return self.getInputFilename(); })
        .def("setTaxonomyFilename", [](OpenMS::XTandemInfile& self, const OpenMS::String& filename) { return self.setTaxonomyFilename(filename); }, "filename"_a)
        .def("getTaxonomyFilename", [](const OpenMS::XTandemInfile& self) { return self.getTaxonomyFilename(); })
        .def("setDefaultParametersFilename", [](OpenMS::XTandemInfile& self, const OpenMS::String& filename) { return self.setDefaultParametersFilename(filename); }, "filename"_a)
        .def("getDefaultParametersFilename", [](const OpenMS::XTandemInfile& self) { return self.getDefaultParametersFilename(); })
        .def("setTaxon", [](OpenMS::XTandemInfile& self, const OpenMS::String& taxon) { return self.setTaxon(taxon); }, "taxon"_a)
        .def("getTaxon", [](const OpenMS::XTandemInfile& self) { return self.getTaxon(); })
        .def("setMaxPrecursorCharge", [](OpenMS::XTandemInfile& self, int max_charge) { return self.setMaxPrecursorCharge(max_charge); }, "max_charge"_a)
        .def("getMaxPrecursorCharge", [](const OpenMS::XTandemInfile& self) { return self.getMaxPrecursorCharge(); })
        .def("setNumberOfMissedCleavages", [](OpenMS::XTandemInfile& self, unsigned int missed_cleavages) { return self.setNumberOfMissedCleavages(missed_cleavages); }, "missed_cleavages"_a)
        .def("getNumberOfMissedCleavages", [](const OpenMS::XTandemInfile& self) { return self.getNumberOfMissedCleavages(); })
        .def("setOutputResults", [](OpenMS::XTandemInfile& self, const OpenMS::String& result) { return self.setOutputResults(result); }, "result"_a)
        .def("getOutputResults", [](const OpenMS::XTandemInfile& self) { return self.getOutputResults(); })
        .def("setMaxValidEValue", [](OpenMS::XTandemInfile& self, double value) { return self.setMaxValidEValue(value); }, "value"_a)
        .def("getMaxValidEValue", [](const OpenMS::XTandemInfile& self) { return self.getMaxValidEValue(); })
        .def("setSemiCleavage", [](OpenMS::XTandemInfile& self, bool semi_cleavage) { return self.setSemiCleavage(semi_cleavage); }, "semi_cleavage"_a)
        .def("setAllowIsotopeError", [](OpenMS::XTandemInfile& self, bool allow_isotope_error) { return self.setAllowIsotopeError(allow_isotope_error); }, "allow_isotope_error"_a)
        .def("setCleavageSite", [](OpenMS::XTandemInfile& self, const OpenMS::String& cleavage_site) { return self.setCleavageSite(cleavage_site); }, "cleavage_site"_a)
        .def("getCleavageSite", [](const OpenMS::XTandemInfile& self) { return self.getCleavageSite(); })
        .def("write", [](OpenMS::XTandemInfile& self, const OpenMS::String& filename, bool ignore_member_parameters, bool force_default_mods) { return self.write(filename, ignore_member_parameters, force_default_mods); }, "filename"_a, "ignore_member_parameters"_a = false, "force_default_mods"_a = false)
        ;
    // ErrorUnit enum nested under XTandemInfile
    nb::enum_<OpenMS::XTandemInfile::ErrorUnit>(xtandeminfile_class, "ErrorUnit", nb::is_arithmetic())
        .value("DALTONS", OpenMS::XTandemInfile::ErrorUnit::DALTONS)
        .value("PPM", OpenMS::XTandemInfile::ErrorUnit::PPM)
        .export_values();
    // MassType enum nested under XTandemInfile
    nb::enum_<OpenMS::XTandemInfile::MassType>(xtandeminfile_class, "MassType", nb::is_arithmetic())
        .value("MONOISOTOPIC", OpenMS::XTandemInfile::MassType::MONOISOTOPIC)
        .value("AVERAGE", OpenMS::XTandemInfile::MassType::AVERAGE)
        .export_values();

    // -----------------------------------------------------------------------
    // XTandemXMLFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::XTandemXMLFile, OpenMS::Internal::XMLFile>(m, "XTandemXMLFile", "Used to load XTandemXML files")
        .def(nb::init<>())
        .def("load", [](OpenMS::XTandemXMLFile& self, const OpenMS::String& filename, OpenMS::ProteinIdentification& protein_identification, OpenMS::PeptideIdentificationList& id_data, OpenMS::ModificationDefinitionsSet& mod_def_set) { return self.load(filename, protein_identification, id_data, mod_def_set); }, "filename"_a, "protein_identification"_a, "id_data"_a, "mod_def_set"_a)
        ;


    // -----------------------------------------------------------------------
    // IsobaricChannelExtractor
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsobaricChannelExtractor, OpenMS::DefaultParamHandler>(m, "IsobaricChannelExtractor", "OpenMS class IsobaricChannelExtractor")
        .def(nb::init<const OpenMS::ItraqFourPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::ItraqEightPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTSixPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTTenPlexQuantitationMethod*>(), "quant_method"_a)
        .def("extractChannels", &OpenMS::IsobaricChannelExtractor::extractChannels, "ms_exp_data"_a, "consensus_map"_a, "Extract isobaric channels")
        ;


    // -----------------------------------------------------------------------
    // DIAScoring
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DIAScoring, OpenMS::DefaultParamHandler>(m, "DIAScoring", "OpenMS class DIAScoring")
        .def(nb::init<>())
        .def("score_with_isotopes", [](OpenMS::DIAScoring& self,
                std::vector<std::shared_ptr<OpenSwath::OSSpectrum>>& spectrum,
                std::vector<OpenSwath::LightTransition>& transitions,
                OpenMS::RangeMobility& im_range,
                double dotprod, double manhattan) {
            self.score_with_isotopes(spectrum, transitions, im_range, dotprod, manhattan);
            return nb::make_tuple(dotprod, manhattan);
        }, "spectrum"_a, "transitions"_a, "im_range"_a, "dotprod"_a, "manhattan"_a)
        .def("dia_by_ion_score", [](OpenMS::DIAScoring& self,
                std::vector<std::shared_ptr<OpenSwath::OSSpectrum>>& spectrum,
                OpenMS::AASequence& sequence,
                int charge,
                OpenMS::RangeMobility& im_range,
                double bseries_score, double yseries_score) {
            self.dia_by_ion_score(spectrum, sequence, charge, im_range, bseries_score, yseries_score);
            return nb::make_tuple(bseries_score, yseries_score);
        }, "spectrum"_a, "sequence"_a, "charge"_a, "im_range"_a, "bseries_score"_a, "yseries_score"_a,
        "Score the DIA window for b/y ion series")
        .def("dia_ms1_massdiff_score", [](const OpenMS::DIAScoring& self, double precursor_mz, const std::vector<OpenSwath::SpectrumPtr>& spectrum, const OpenMS::RangeMobility& im_range) { double ppm_score; bool result = self.dia_ms1_massdiff_score(precursor_mz, spectrum, im_range, ppm_score); return nb::make_tuple(result, ppm_score); }, "precursor_mz"_a, "spectrum"_a, "im_range"_a, "Score MS1 mass difference")
        .def("dia_ms1_isotope_scores_averagine", [](const OpenMS::DIAScoring& self,
                double precursor_mz,
                const std::vector<OpenSwath::SpectrumPtr>& spectrum,
                int charge_state,
                OpenMS::RangeMobility& im_range) {
            double isotope_corr = 0, isotope_overlap = 0;
            self.dia_ms1_isotope_scores_averagine(precursor_mz, spectrum, charge_state, im_range, isotope_corr, isotope_overlap);
            return nb::make_tuple(isotope_corr, isotope_overlap);
        }, "precursor_mz"_a, "spectrum"_a, "charge_state"_a, "im_range"_a,
        "Precursor isotope scores using averagine model. Returns (isotope_corr, isotope_overlap)")
        .def("dia_ms1_isotope_scores", [](const OpenMS::DIAScoring& self,
                double precursor_mz,
                const std::vector<OpenSwath::SpectrumPtr>& spectrum,
                OpenMS::RangeMobility& im_range,
                const OpenMS::EmpiricalFormula& sum_formula) {
            double isotope_corr = 0, isotope_overlap = 0;
            self.dia_ms1_isotope_scores(precursor_mz, spectrum, im_range, isotope_corr, isotope_overlap, sum_formula);
            return nb::make_tuple(isotope_corr, isotope_overlap);
        }, "precursor_mz"_a, "spectrum"_a, "im_range"_a, "sum_formula"_a,
        "Precursor isotope scores using empirical formula. Returns (isotope_corr, isotope_overlap)")
        ;


    // -----------------------------------------------------------------------
    // MRMFeatureFinderScoring
    // -----------------------------------------------------------------------
    auto mrmfeaturefinderscoring_class = nb::class_<OpenMS::MRMFeatureFinderScoring, OpenMS::DefaultParamHandler>(m, "MRMFeatureFinderScoring", "OpenMS class MRMFeatureFinderScoring")
        .def(nb::init<>())
        .def("pickExperiment", [](OpenMS::MRMFeatureFinderScoring& self,
                OpenMS::MSExperiment& chromatograms,
                OpenMS::FeatureMap& output,
                OpenMS::TargetedExperiment& transition_exp,
                OpenMS::TransformationDescription& trafo,
                OpenMS::MSExperiment& swath_map) {
            self.pickExperiment(chromatograms, output, transition_exp, trafo, swath_map);
        }, "chromatograms"_a, "output"_a, "transition_exp"_a, "trafo"_a, "swath_map"_a)
        .def("setStrictFlag", &OpenMS::MRMFeatureFinderScoring::setStrictFlag, "flag"_a)
        .def("setMS1Map", [](OpenMS::MRMFeatureFinderScoring& self, std::shared_ptr<OpenMS::SpectrumAccessOpenMS> ms1_map) {
            self.setMS1Map(ms1_map);
        }, "ms1_map"_a)
        .def("prepareProteinPeptideMaps_", [](OpenMS::MRMFeatureFinderScoring& self, const OpenSwath::LightTargetedExperiment& transition_exp) { self.prepareProteinPeptideMaps_(transition_exp); }, "transition_exp"_a, "Prepares the internal mappings of peptides and proteins")
        .def("scorePeakgroups", [](const OpenMS::MRMFeatureFinderScoring& self,
                OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& transition_group,
                const OpenMS::TransformationDescription& trafo,
                const std::vector<OpenSwath::SwathMap>& swath_maps,
                OpenMS::FeatureMap& output,
                bool ms1only) {
            self.scorePeakgroups(transition_group, trafo, swath_maps, output, ms1only);
        }, "transition_group"_a, "trafo"_a, "swath_maps"_a, "output"_a, "ms1only"_a = false,
           "Score peak groups in a transition group")
        ;
    def_ProgressLogger<OpenMS::MRMFeatureFinderScoring>(mrmfeaturefinderscoring_class);

    // Free function aliases for backward compatibility
    m.def("isPEFFFile", [](const OpenMS::String& filename) {
        return OpenMS::PEFFFile::isPEFFFile(filename);
    }, "filename"_a, "Check if a file is in PEFF format");
    m.def("toProForma", [](const OpenMS::PEFFEntry& entry) {
        return OpenMS::PEFFFile::toProForma(entry);
    }, "entry"_a, "Convert a PEFFEntry to ProForma notation");

    // -----------------------------------------------------------------------
    // __static_* module-level wrappers for IMTypes
    // -----------------------------------------------------------------------
    m.def("__static_IMTypes_determineIMFormat", [](const OpenMS::MSExperiment& exp, int ms_level) -> OpenMS::IMFormat { return OpenMS::IMTypes::determineIMFormat(exp, ms_level); }, "exp"_a, "ms_level"_a);
    m.def("__static_IMTypes_toDriftTimeUnit", [](const OpenMS::String& dtu_string) -> OpenMS::DriftTimeUnit { return OpenMS::toDriftTimeUnit(dtu_string); }, "dtu_string"_a);
    m.def("__static_IMTypes_driftTimeUnitToString", [](OpenMS::DriftTimeUnit value) -> OpenMS::String { return OpenMS::driftTimeUnitToString(value); }, "value"_a);
    m.def("__static_IMTypes_toIMFormat", [](const OpenMS::String& im_format) -> OpenMS::IMFormat { return OpenMS::toIMFormat(im_format); }, "im_format"_a);
    m.def("__static_IMTypes_toIMPeakType", [](const OpenMS::String& im_peak_type) -> OpenMS::IMPeakType { return OpenMS::toIMPeakType(im_peak_type); }, "im_peak_type"_a);
    m.def("__static_IMTypes_imPeakTypeToString", [](OpenMS::IMPeakType value) -> OpenMS::String { return OpenMS::imPeakTypeToString(value); }, "value"_a);

}
