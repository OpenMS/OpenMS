// pyOpenMS nanobind bindings
// Domain: analysis

#include "all_casters.h"
#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/ANALYSIS/ID/FIAMSScheduler.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>
#include <OpenMS/ANALYSIS/ID/MorpheusScore.h>
#include <OpenMS/ANALYSIS/ID/PScore.h>
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
// #include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>
// #include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeaturePicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>
#include <OpenMS/ANALYSIS/XLMS/OpenPepXLAlgorithm.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/ANALYSIS/XLMS/XQuestScores.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
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

NB_MODULE(_pyopenms_analysis, m) {
    m.doc() = "pyOpenMS analysis bindings";

    // -----------------------------------------------------------------------
    // OpenSwath_Scores
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSwath_Scores>(m, "OpenSwath_Scores", "OpenSWATH scoring results struct")
        .def(nb::init<>())
        .def_rw("elution_model_fit_score", &OpenMS::OpenSwath_Scores::elution_model_fit_score)
        .def_rw("library_corr", &OpenMS::OpenSwath_Scores::library_corr)
        .def_rw("library_norm_manhattan", &OpenMS::OpenSwath_Scores::library_norm_manhattan)
        .def_rw("library_rootmeansquare", &OpenMS::OpenSwath_Scores::library_rootmeansquare)
        .def_rw("library_sangle", &OpenMS::OpenSwath_Scores::library_sangle)
        .def_rw("norm_rt_score", &OpenMS::OpenSwath_Scores::norm_rt_score)
        .def_rw("isotope_correlation", &OpenMS::OpenSwath_Scores::isotope_correlation)
        .def_rw("isotope_overlap", &OpenMS::OpenSwath_Scores::isotope_overlap)
        .def_rw("massdev_score", &OpenMS::OpenSwath_Scores::massdev_score)
        .def_rw("xcorr_coelution_score", &OpenMS::OpenSwath_Scores::xcorr_coelution_score)
        .def_rw("xcorr_shape_score", &OpenMS::OpenSwath_Scores::xcorr_shape_score)
        .def_rw("yseries_score", &OpenMS::OpenSwath_Scores::yseries_score)
        .def_rw("bseries_score", &OpenMS::OpenSwath_Scores::bseries_score)
        .def_rw("log_sn_score", &OpenMS::OpenSwath_Scores::log_sn_score)
        .def_rw("weighted_coelution_score", &OpenMS::OpenSwath_Scores::weighted_coelution_score)
        .def_rw("weighted_xcorr_shape", &OpenMS::OpenSwath_Scores::weighted_xcorr_shape)
        .def_rw("weighted_massdev_score", &OpenMS::OpenSwath_Scores::weighted_massdev_score)
        .def_rw("intensity", &OpenMS::OpenSwath_Scores::intensity)
        .def_rw("total_xic", &OpenMS::OpenSwath_Scores::total_xic)
        .def_rw("nr_peaks", &OpenMS::OpenSwath_Scores::nr_peaks)
        .def_rw("sn_ratio", &OpenMS::OpenSwath_Scores::sn_ratio)
        .def_rw("mi_score", &OpenMS::OpenSwath_Scores::mi_score)
        .def_rw("weighted_mi_score", &OpenMS::OpenSwath_Scores::weighted_mi_score)
        .def_rw("rt_difference", &OpenMS::OpenSwath_Scores::rt_difference)
        .def_rw("normalized_experimental_rt", &OpenMS::OpenSwath_Scores::normalized_experimental_rt)
        .def_rw("raw_rt_score", &OpenMS::OpenSwath_Scores::raw_rt_score)
        // ms1 scores
        .def_rw("ms1_xcorr_coelution_score", &OpenMS::OpenSwath_Scores::ms1_xcorr_coelution_score)
        .def_rw("ms1_xcorr_coelution_contrast_score", &OpenMS::OpenSwath_Scores::ms1_xcorr_coelution_contrast_score)
        .def_rw("ms1_xcorr_coelution_combined_score", &OpenMS::OpenSwath_Scores::ms1_xcorr_coelution_combined_score)
        .def_rw("ms1_xcorr_shape_score", &OpenMS::OpenSwath_Scores::ms1_xcorr_shape_score)
        .def_rw("ms1_xcorr_shape_contrast_score", &OpenMS::OpenSwath_Scores::ms1_xcorr_shape_contrast_score)
        .def_rw("ms1_xcorr_shape_combined_score", &OpenMS::OpenSwath_Scores::ms1_xcorr_shape_combined_score)
        .def_rw("ms1_ppm_score", &OpenMS::OpenSwath_Scores::ms1_ppm_score)
        .def_rw("ms1_isotope_correlation", &OpenMS::OpenSwath_Scores::ms1_isotope_correlation)
        .def_rw("ms1_isotope_overlap", &OpenMS::OpenSwath_Scores::ms1_isotope_overlap)
        .def_rw("ms1_mi_score", &OpenMS::OpenSwath_Scores::ms1_mi_score)
        .def_rw("ms1_mi_contrast_score", &OpenMS::OpenSwath_Scores::ms1_mi_contrast_score)
        .def_rw("ms1_mi_combined_score", &OpenMS::OpenSwath_Scores::ms1_mi_combined_score)
        // IM scores
        .def_rw("im_xcorr_coelution_score", &OpenMS::OpenSwath_Scores::im_xcorr_coelution_score)
        .def_rw("im_xcorr_shape_score", &OpenMS::OpenSwath_Scores::im_xcorr_shape_score)
        .def_rw("im_delta_score", &OpenMS::OpenSwath_Scores::im_delta_score)
        .def_rw("im_ms1_delta_score", &OpenMS::OpenSwath_Scores::im_ms1_delta_score)
        .def_rw("im_drift", &OpenMS::OpenSwath_Scores::im_drift)
        .def_rw("im_drift_left", &OpenMS::OpenSwath_Scores::im_drift_left)
        .def_rw("im_drift_right", &OpenMS::OpenSwath_Scores::im_drift_right)
        .def_rw("im_drift_weighted", &OpenMS::OpenSwath_Scores::im_drift_weighted)
        .def_rw("im_delta", &OpenMS::OpenSwath_Scores::im_delta)
        .def_rw("im_log_intensity", &OpenMS::OpenSwath_Scores::im_log_intensity)
        .def_rw("im_ms1_contrast_coelution", &OpenMS::OpenSwath_Scores::im_ms1_contrast_coelution)
        .def_rw("im_ms1_contrast_shape", &OpenMS::OpenSwath_Scores::im_ms1_contrast_shape)
        .def_rw("im_ms1_sum_contrast_coelution", &OpenMS::OpenSwath_Scores::im_ms1_sum_contrast_coelution)
        .def_rw("im_ms1_sum_contrast_shape", &OpenMS::OpenSwath_Scores::im_ms1_sum_contrast_shape)
        .def_rw("im_ms1_drift", &OpenMS::OpenSwath_Scores::im_ms1_drift)
        .def_rw("im_ms1_delta", &OpenMS::OpenSwath_Scores::im_ms1_delta)
        .def_rw("im_ind_contrast_coelution", &OpenMS::OpenSwath_Scores::im_ind_contrast_coelution)
        .def_rw("im_ind_contrast_shape", &OpenMS::OpenSwath_Scores::im_ind_contrast_shape)
        .def_rw("im_ind_sum_contrast_coelution", &OpenMS::OpenSwath_Scores::im_ind_sum_contrast_coelution)
        .def_rw("im_ind_sum_contrast_shape", &OpenMS::OpenSwath_Scores::im_ind_sum_contrast_shape)
        // DIA scores
        .def_rw("dotprod_score_dia", &OpenMS::OpenSwath_Scores::dotprod_score_dia)
        .def_rw("manhatt_score_dia", &OpenMS::OpenSwath_Scores::manhatt_score_dia)
        .def_rw("library_manhattan", &OpenMS::OpenSwath_Scores::library_manhattan)
        .def_rw("library_dotprod", &OpenMS::OpenSwath_Scores::library_dotprod)
        .def("calculate_lda_prescore", [](const OpenMS::OpenSwath_Scores& self, const OpenMS::OpenSwath_Scores& scores) { return self.calculate_lda_prescore(scores); }, "scores"_a, "Calculate LDA prescore")
        .def("calculate_swath_lda_prescore", [](const OpenMS::OpenSwath_Scores& self, const OpenMS::OpenSwath_Scores& scores) { return self.calculate_swath_lda_prescore(scores); }, "scores"_a, "Calculate SWATH LDA prescore")
        .def("get_quick_lda_score", [](const OpenMS::OpenSwath_Scores& self, double library_corr_, double library_norm_manhattan_, double norm_rt_score_, double xcorr_coelution_score_, double xcorr_shape_score_, double log_sn_score_) { return self.get_quick_lda_score(library_corr_, library_norm_manhattan_, norm_rt_score_, xcorr_coelution_score_, xcorr_shape_score_, log_sn_score_); }, "library_corr"_a, "library_norm_manhattan"_a, "norm_rt_score"_a, "xcorr_coelution_score"_a, "xcorr_shape_score"_a, "log_sn_score"_a, "Get quick LDA score")
        ;

    // -----------------------------------------------------------------------
    // OpenSwath_Scores_Usage
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSwath_Scores_Usage>(m, "OpenSwath_Scores_Usage", "OpenSWATH scores usage flags")
        .def(nb::init<>())
        .def_rw("use_coelution_score_", &OpenMS::OpenSwath_Scores_Usage::use_coelution_score_)
        .def_rw("use_shape_score_", &OpenMS::OpenSwath_Scores_Usage::use_shape_score_)
        .def_rw("use_rt_score_", &OpenMS::OpenSwath_Scores_Usage::use_rt_score_)
        .def_rw("use_library_score_", &OpenMS::OpenSwath_Scores_Usage::use_library_score_)
        .def_rw("use_elution_model_score_", &OpenMS::OpenSwath_Scores_Usage::use_elution_model_score_)
        .def_rw("use_intensity_score_", &OpenMS::OpenSwath_Scores_Usage::use_intensity_score_)
        .def_rw("use_total_xic_score_", &OpenMS::OpenSwath_Scores_Usage::use_total_xic_score_)
        .def_rw("use_total_mi_score_", &OpenMS::OpenSwath_Scores_Usage::use_total_mi_score_)
        .def_rw("use_nr_peaks_score_", &OpenMS::OpenSwath_Scores_Usage::use_nr_peaks_score_)
        .def_rw("use_sn_score_", &OpenMS::OpenSwath_Scores_Usage::use_sn_score_)
        .def_rw("use_mi_score_", &OpenMS::OpenSwath_Scores_Usage::use_mi_score_)
        .def_rw("use_dia_scores_", &OpenMS::OpenSwath_Scores_Usage::use_dia_scores_)
        .def_rw("use_im_scores", &OpenMS::OpenSwath_Scores_Usage::use_im_scores)
        .def_rw("use_ms1_correlation", &OpenMS::OpenSwath_Scores_Usage::use_ms1_correlation)
        .def_rw("use_ms1_fullscan", &OpenMS::OpenSwath_Scores_Usage::use_ms1_fullscan)
        .def_rw("use_ms1_mi", &OpenMS::OpenSwath_Scores_Usage::use_ms1_mi)
        .def_rw("use_uis_scores", &OpenMS::OpenSwath_Scores_Usage::use_uis_scores)
        .def_rw("use_ionseries_scores", &OpenMS::OpenSwath_Scores_Usage::use_ionseries_scores)
        .def_rw("use_ms2_isotope_scores", &OpenMS::OpenSwath_Scores_Usage::use_ms2_isotope_scores)
        .def_rw("use_peak_shape_metrics", &OpenMS::OpenSwath_Scores_Usage::use_peak_shape_metrics)
        ;

    // -----------------------------------------------------------------------
    // AbsoluteQuantitationMethod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitationMethod>(m, "AbsoluteQuantitationMethod", 
        R"doc(
Holds information about a quantitation method for absolute quantitation using
Isotope Dilution Mass Spectrometry (IDMS). This includes calibration curve
parameters, limits of detection (LOD), limits of quantitation (LOQ), and
the transformation model used for concentration calculation
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::AbsoluteQuantitationMethod &>())
        .def("__copy__", [](const OpenMS::AbsoluteQuantitationMethod& self) { return OpenMS::AbsoluteQuantitationMethod(self); })
        .def("__deepcopy__", [](const OpenMS::AbsoluteQuantitationMethod& self, nb::dict) { return OpenMS::AbsoluteQuantitationMethod(self); }, "memo"_a)
        .def("setComponentName", [](OpenMS::AbsoluteQuantitationMethod& self, const OpenMS::String& component_name) { return self.setComponentName(component_name); }, "component_name"_a, "Sets the component name (unique identifier for the analyte)")
        .def("getComponentName", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getComponentName(); }, "Returns the component name (unique identifier for the analyte)")
        .def("setFeatureName", [](OpenMS::AbsoluteQuantitationMethod& self, const OpenMS::String& feature_name) { return self.setFeatureName(feature_name); }, "feature_name"_a, "Sets the feature name (e.g., peak_apex_int or peak_area)")
        .def("getFeatureName", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getFeatureName(); }, "Returns the feature name (e.g., peak_apex_int or peak_area)")
        .def("setISName", [](OpenMS::AbsoluteQuantitationMethod& self, const OpenMS::String& IS_name) { return self.setISName(IS_name); }, "IS_name"_a, "Sets the internal standard name used for ratio calculation")
        .def("getISName", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getISName(); }, "Returns the internal standard name used for ratio calculation")
        .def("setLLOD", [](OpenMS::AbsoluteQuantitationMethod& self, double llod) { return self.setLLOD(llod); }, "llod"_a, "Sets the lower limit of detection (LLOD)")
        .def("getLLOD", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getLLOD(); }, "Returns the lower limit of detection (LLOD)")
        .def("setULOD", [](OpenMS::AbsoluteQuantitationMethod& self, double ulod) { return self.setULOD(ulod); }, "ulod"_a, "Sets the upper limit of detection (ULOD)")
        .def("getULOD", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getULOD(); }, "Returns the upper limit of detection (ULOD)")
        .def("checkLOD", [](const OpenMS::AbsoluteQuantitationMethod& self, double value) { return self.checkLOD(value); }, "value"_a, "Checks if the given value is within the limits of detection (LOD)")
        .def("setLLOQ", [](OpenMS::AbsoluteQuantitationMethod& self, double lloq) { return self.setLLOQ(lloq); }, "lloq"_a, "Sets the lower limit of quantitation (LLOQ)")
        .def("getLLOQ", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getLLOQ(); }, "Returns the lower limit of quantitation (LLOQ)")
        .def("setULOQ", [](OpenMS::AbsoluteQuantitationMethod& self, double uloq) { return self.setULOQ(uloq); }, "uloq"_a, "Sets the upper limit of quantitation (ULOQ)")
        .def("getULOQ", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getULOQ(); }, "Returns the upper limit of quantitation (ULOQ)")
        .def("checkLOQ", [](const OpenMS::AbsoluteQuantitationMethod& self, double value) { return self.checkLOQ(value); }, "value"_a, "Checks if the given value is within the limits of quantitation (LOQ)")
        .def("setNPoints", [](OpenMS::AbsoluteQuantitationMethod& self, int n_points) { return self.setNPoints(n_points); }, "n_points"_a, "Sets the number of points used in the calibration curve")
        .def("getNPoints", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getNPoints(); }, "Returns the number of points used in the calibration curve")
        .def("setCorrelationCoefficient", [](OpenMS::AbsoluteQuantitationMethod& self, double correlation_coefficient) { return self.setCorrelationCoefficient(correlation_coefficient); }, "correlation_coefficient"_a, "Sets the Pearson correlation coefficient of the calibration curve fit")
        .def("getCorrelationCoefficient", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getCorrelationCoefficient(); }, "Returns the Pearson correlation coefficient of the calibration curve fit")
        .def("setConcentrationUnits", [](OpenMS::AbsoluteQuantitationMethod& self, const OpenMS::String& concentration_units) { return self.setConcentrationUnits(concentration_units); }, "concentration_units"_a, "Sets the concentration units for the component")
        .def("getConcentrationUnits", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getConcentrationUnits(); }, "Returns the concentration units for the component")
        .def("setTransformationModel", [](OpenMS::AbsoluteQuantitationMethod& self, const OpenMS::String& transformation_model) { return self.setTransformationModel(transformation_model); }, "transformation_model"_a, "Sets the transformation model type (e.g., linear, b_spline)")
        .def("getTransformationModel", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getTransformationModel(); }, "Returns the transformation model type (e.g., linear, b_spline)")
        .def("setTransformationModelParams", [](OpenMS::AbsoluteQuantitationMethod& self, const OpenMS::Param& transformation_model_params) { return self.setTransformationModelParams(transformation_model_params); }, "transformation_model_params"_a, "Sets the transformation model parameters")
        .def("getTransformationModelParams", [](const OpenMS::AbsoluteQuantitationMethod& self) { return self.getTransformationModelParams(); }, "Returns the transformation model parameters")
        ;

    // -----------------------------------------------------------------------
    // AccurateMassSearchResult
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AccurateMassSearchResult>(m, "AccurateMassSearchResult", "OpenMS class AccurateMassSearchResult")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::AccurateMassSearchResult &>())
        .def("__copy__", [](const OpenMS::AccurateMassSearchResult& self) { return OpenMS::AccurateMassSearchResult(self); })
        .def("__deepcopy__", [](const OpenMS::AccurateMassSearchResult& self, nb::dict) { return OpenMS::AccurateMassSearchResult(self); }, "memo"_a)
        .def("getObservedMZ", [](const OpenMS::AccurateMassSearchResult& self) { return self.getObservedMZ(); })
        .def("setObservedMZ", [](OpenMS::AccurateMassSearchResult& self, const double& p0) { return self.setObservedMZ(p0); })
        .def("getCalculatedMZ", [](const OpenMS::AccurateMassSearchResult& self) { return self.getCalculatedMZ(); })
        .def("setCalculatedMZ", [](OpenMS::AccurateMassSearchResult& self, const double& p0) { return self.setCalculatedMZ(p0); })
        .def("getQueryMass", [](const OpenMS::AccurateMassSearchResult& self) { return self.getQueryMass(); })
        .def("setQueryMass", [](OpenMS::AccurateMassSearchResult& self, const double& p0) { return self.setQueryMass(p0); })
        .def("getFoundMass", [](const OpenMS::AccurateMassSearchResult& self) { return self.getFoundMass(); })
        .def("setFoundMass", [](OpenMS::AccurateMassSearchResult& self, const double& p0) { return self.setFoundMass(p0); })
        .def("getCharge", [](const OpenMS::AccurateMassSearchResult& self) { return self.getCharge(); })
        .def("setCharge", [](OpenMS::AccurateMassSearchResult& self, const int& p0) { return self.setCharge(p0); })
        .def("getMZErrorPPM", [](const OpenMS::AccurateMassSearchResult& self) { return self.getMZErrorPPM(); })
        .def("setMZErrorPPM", [](OpenMS::AccurateMassSearchResult& self, double p0) { return self.setMZErrorPPM(p0); })
        .def("getObservedRT", [](const OpenMS::AccurateMassSearchResult& self) { return self.getObservedRT(); })
        .def("setObservedRT", [](OpenMS::AccurateMassSearchResult& self, const double& rt) { return self.setObservedRT(rt); }, "rt"_a)
        .def("getObservedIntensity", [](const OpenMS::AccurateMassSearchResult& self) { return self.getObservedIntensity(); })
        .def("setObservedIntensity", [](OpenMS::AccurateMassSearchResult& self, const double& p0) { return self.setObservedIntensity(p0); })
        .def("getIndividualIntensities", [](const OpenMS::AccurateMassSearchResult& self) { return self.getIndividualIntensities(); })
        .def("getMatchingIndex", [](const OpenMS::AccurateMassSearchResult& self) { return self.getMatchingIndex(); })
        .def("setMatchingIndex", [](OpenMS::AccurateMassSearchResult& self, const size_t& p0) { return self.setMatchingIndex(p0); })
        .def("getSourceFeatureIndex", [](const OpenMS::AccurateMassSearchResult& self) { return self.getSourceFeatureIndex(); })
        .def("setSourceFeatureIndex", [](OpenMS::AccurateMassSearchResult& self, const size_t& p0) { return self.setSourceFeatureIndex(p0); })
        .def("getFoundAdduct", [](const OpenMS::AccurateMassSearchResult& self) { return self.getFoundAdduct(); })
        .def("setFoundAdduct", [](OpenMS::AccurateMassSearchResult& self, const OpenMS::String& p0) { return self.setFoundAdduct(p0); })
        .def("getFormulaString", [](const OpenMS::AccurateMassSearchResult& self) { return self.getFormulaString(); })
        .def("setEmpiricalFormula", [](OpenMS::AccurateMassSearchResult& self, const OpenMS::String& p0) { return self.setEmpiricalFormula(p0); })
        .def("getMatchingHMDBids", [](const OpenMS::AccurateMassSearchResult& self) -> const std::vector<OpenMS::String> & { return self.getMatchingHMDBids(); }, nb::rv_policy::reference_internal)
        .def("setMatchingHMDBids", [](OpenMS::AccurateMassSearchResult& self, const std::vector<OpenMS::String>& p0) { return self.setMatchingHMDBids(p0); })
        .def("getMasstraceIntensities", [](const OpenMS::AccurateMassSearchResult& self) -> const std::vector<double> & { return self.getMasstraceIntensities(); }, nb::rv_policy::reference_internal)
        .def("setMasstraceIntensities", [](OpenMS::AccurateMassSearchResult& self, const std::vector<double>& p0) { return self.setMasstraceIntensities(p0); })
        .def("getIsotopesSimScore", [](const OpenMS::AccurateMassSearchResult& self) { return self.getIsotopesSimScore(); })
        .def("setIsotopesSimScore", [](OpenMS::AccurateMassSearchResult& self, const double& p0) { return self.setIsotopesSimScore(p0); })
        .def("setIndividualIntensities", [](OpenMS::AccurateMassSearchResult& self, const std::vector<double>& intensities) { self.setIndividualIntensities(intensities); }, "intensities"_a, "Sets the individual intensities")
        ;

    // -----------------------------------------------------------------------
    // CV
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::CV>(m, "CV", "OpenMS class CV")
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::String>())
        .def("__copy__", [](const OpenMS::TargetedExperimentHelper::CV& self) { return OpenMS::TargetedExperimentHelper::CV(self); })
        .def("__deepcopy__", [](const OpenMS::TargetedExperimentHelper::CV& self, nb::dict) { return OpenMS::TargetedExperimentHelper::CV(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def_rw("id", &OpenMS::TargetedExperimentHelper::CV::id)
        .def_rw("fullname", &OpenMS::TargetedExperimentHelper::CV::fullname)
        .def_rw("version", &OpenMS::TargetedExperimentHelper::CV::version)
        .def_rw("URI", &OpenMS::TargetedExperimentHelper::CV::URI)
        ;

    // -----------------------------------------------------------------------
    // ClusterProxyKD
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ClusterProxyKD>(m, "ClusterProxyKD", "Proxy for a (potential) cluster")
        .def(nb::init<>())
        .def(nb::init<size_t, double, size_t>())
        .def(nb::init<const OpenMS::ClusterProxyKD &>())
        .def("__copy__", [](const OpenMS::ClusterProxyKD& self) { return OpenMS::ClusterProxyKD(self); })
        .def("__deepcopy__", [](const OpenMS::ClusterProxyKD& self, nb::dict) { return OpenMS::ClusterProxyKD(self); }, "memo"_a)
        .def(nb::self < nb::self)
        .def(nb::self != nb::self)
        .def(nb::self == nb::self)
        .def("getSize", [](const OpenMS::ClusterProxyKD& self) { return self.getSize(); })
        .def("isValid", [](const OpenMS::ClusterProxyKD& self) { return self.isValid(); })
        .def("getAvgDistance", [](const OpenMS::ClusterProxyKD& self) { return self.getAvgDistance(); })
        .def("getCenterIndex", [](const OpenMS::ClusterProxyKD& self) { return self.getCenterIndex(); })
        ;

    // -----------------------------------------------------------------------
    // Compound
    // -----------------------------------------------------------------------
    auto compound_class = nb::class_<OpenMS::TargetedExperimentHelper::Compound>(m, "Compound", 
        R"doc(
Represents a compound in a targeted experiment (e.g. used by
ReactionMonitoringTransition and IncludeExcludeTarget)
CVTermList
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TargetedExperimentHelper::Compound &>())
        .def("__copy__", [](const OpenMS::TargetedExperimentHelper::Compound& self) { return OpenMS::TargetedExperimentHelper::Compound(self); })
        .def("__deepcopy__", [](const OpenMS::TargetedExperimentHelper::Compound& self, nb::dict) { return OpenMS::TargetedExperimentHelper::Compound(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("setChargeState", [](OpenMS::TargetedExperimentHelper::Compound& self, int charge) { return self.setChargeState(charge); }, "charge"_a, "Sets the peptide or compound charge state")
        .def("hasCharge", [](const OpenMS::TargetedExperimentHelper::Compound& self) { return self.hasCharge(); }, "Whether peptide or compound has set charge state")
        .def("getChargeState", [](const OpenMS::TargetedExperimentHelper::Compound& self) { return self.getChargeState(); }, "Returns the peptide or compound charge state")
        .def("hasRetentionTime", [](const OpenMS::TargetedExperimentHelper::Compound& self) { return self.hasRetentionTime(); }, "Check whether compound or peptide has an annotated retention time")
        .def("getRetentionTime", [](const OpenMS::TargetedExperimentHelper::Compound& self) { return self.getRetentionTime(); }, "Gets compound or peptide retention time")
        .def("getRetentionTimeType", [](const OpenMS::TargetedExperimentHelper::Compound& self) { return self.getRetentionTimeType(); }, "Get compound or peptide retentiontime type")
        .def("getRetentionTimeUnit", [](const OpenMS::TargetedExperimentHelper::Compound& self) { return self.getRetentionTimeUnit(); }, "Get compound or peptide retentiontime type")
        .def(nb::self != nb::self)
        
        .def_rw("molecular_formula", &OpenMS::TargetedExperimentHelper::Compound::molecular_formula)
        .def_rw("smiles_string", &OpenMS::TargetedExperimentHelper::Compound::smiles_string)
        .def_rw("theoretical_mass", &OpenMS::TargetedExperimentHelper::Compound::theoretical_mass)
        .def_rw("rts", &OpenMS::TargetedExperimentHelper::Compound::rts)
        .def_rw("id", &OpenMS::TargetedExperimentHelper::Compound::id)
        ;
    def_CVTermList<OpenMS::TargetedExperimentHelper::Compound>(compound_class);
    def_MetaInfoInterface<OpenMS::TargetedExperimentHelper::Compound>(compound_class);

    // -----------------------------------------------------------------------
    // ConsensusMapNormalizerAlgorithmMedian
    // -----------------------------------------------------------------------
    auto consensusmapnormalizeralgorithmmedian_class = nb::class_<OpenMS::ConsensusMapNormalizerAlgorithmMedian>(m, "ConsensusMapNormalizerAlgorithmMedian", "Algorithms of ConsensusMapNormalizer *")
        .def(nb::init<>())
        .def_static("normalizeMaps", [](OpenMS::ConsensusMap& map, OpenMS::ConsensusMapNormalizerAlgorithmMedian::NormalizationMethod method, const OpenMS::String& acc_filter, const OpenMS::String& desc_filter) { return OpenMS::ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, method, acc_filter, desc_filter); }, "map"_a, "method"_a, "acc_filter"_a, "desc_filter"_a, "Normalizes the maps of the consensusMap")
        .def_static("computeMedians", [](const OpenMS::ConsensusMap& map, const OpenMS::String& acc_filter, const OpenMS::String& desc_filter) {
            std::vector<double> medians;
            auto idx = OpenMS::ConsensusMapNormalizerAlgorithmMedian::computeMedians(map, medians, acc_filter, desc_filter);
            return nb::make_tuple(idx, medians);
        }, "map"_a, "acc_filter"_a, "desc_filter"_a, "Computes medians of all maps and returns tuple of (index of map with most features, medians vector)")
        ;
    // NormalizationMethod enum nested under ConsensusMapNormalizerAlgorithmMedian
    nb::enum_<OpenMS::ConsensusMapNormalizerAlgorithmMedian::NormalizationMethod>(consensusmapnormalizeralgorithmmedian_class, "NormalizationMethod", nb::is_arithmetic(), "Normalization method for consensus map alignment (scale or shift)")
        .value("NM_SCALE", OpenMS::ConsensusMapNormalizerAlgorithmMedian::NormalizationMethod::NM_SCALE)
        .value("NM_SHIFT", OpenMS::ConsensusMapNormalizerAlgorithmMedian::NormalizationMethod::NM_SHIFT)

        .export_values();

    // -----------------------------------------------------------------------
    // ConsensusMapNormalizerAlgorithmQuantile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusMapNormalizerAlgorithmQuantile>(m, "ConsensusMapNormalizerAlgorithmQuantile", "Algorithms of ConsensusMapNormalizer *")
        .def(nb::init<>())
        .def_static("normalizeMaps", [](OpenMS::ConsensusMap& map) { return OpenMS::ConsensusMapNormalizerAlgorithmQuantile::normalizeMaps(map); }, "map"_a)
        .def_static("resample", [](const std::vector<double>& data_in, unsigned int n_resampling_points) {
            std::vector<double> data_out;
            OpenMS::ConsensusMapNormalizerAlgorithmQuantile::resample(data_in, data_out, n_resampling_points);
            return data_out;
        }, "data_in"_a, "n_resampling_points"_a, "Resamples data_in and returns the resampled data")
        .def_static("extractIntensityVectors", [](const OpenMS::ConsensusMap& map) { std::vector<std::vector<double>> out_intensities; OpenMS::ConsensusMapNormalizerAlgorithmQuantile::extractIntensityVectors(map, out_intensities); return out_intensities; }, "map"_a, "Extracts the intensities of the features of the different maps")
        .def_static("setNormalizedIntensityValues", [](const std::vector<std::vector<double>>& feature_ints, OpenMS::ConsensusMap& map) { return OpenMS::ConsensusMapNormalizerAlgorithmQuantile::setNormalizedIntensityValues(feature_ints, map); }, "feature_ints"_a, "map"_a, "Sets the normalized intensity values to the consensus map")
        ;

    // -----------------------------------------------------------------------
    // ConsensusMapNormalizerAlgorithmThreshold
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusMapNormalizerAlgorithmThreshold>(m, "ConsensusMapNormalizerAlgorithmThreshold", "Algorithms of ConsensusMapNormalizer *")
        .def(nb::init<>())
        .def_static("computeCorrelation", [](const OpenMS::ConsensusMap& map, const double& ratio_threshold, const OpenMS::String& acc_filter, const OpenMS::String& desc_filter) { return OpenMS::ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(map, ratio_threshold, acc_filter, desc_filter); }, "map"_a, "ratio_threshold"_a, "acc_filter"_a, "desc_filter"_a, "Determines the ratio of all maps to the map with the most features")
        .def_static("normalizeMaps", [](OpenMS::ConsensusMap& map, const std::vector<double>& ratios) { return OpenMS::ConsensusMapNormalizerAlgorithmThreshold::normalizeMaps(map, ratios); }, "map"_a, "ratios"_a, "Applies the given ratio to the maps of the consensusMap")
        ;

    // -----------------------------------------------------------------------
    // DeconvolvedSpectrum
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DeconvolvedSpectrum>(m, "DeconvolvedSpectrum", 
        R"doc(
A class representing a deconvolved spectrum.
DeconvolvedSpectrum consists of PeakGroup instances representing masses.
For MSn n>1, a PeakGroup representing the precursor mass is also added.
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<int>())
        .def(nb::init<const OpenMS::DeconvolvedSpectrum &>())
        .def("__copy__", [](const OpenMS::DeconvolvedSpectrum& self) { return OpenMS::DeconvolvedSpectrum(self); })
        .def("__deepcopy__", [](const OpenMS::DeconvolvedSpectrum& self, nb::dict) { return OpenMS::DeconvolvedSpectrum(self); }, "memo"_a)
        .def("toSpectrum", [](OpenMS::DeconvolvedSpectrum& self, int to_charge, double tol, bool retain_undeconvolved) { return self.toSpectrum(to_charge, tol, retain_undeconvolved); }, "to_charge"_a, "tol"_a = 10.0, "retain_undeconvolved"_a = false)
        .def("getOriginalSpectrum", [](const OpenMS::DeconvolvedSpectrum& self) -> const OpenMS::MSSpectrum & { return self.getOriginalSpectrum(); }, nb::rv_policy::reference_internal, "Returns the original spectrum")
        .def("getPrecursorPeakGroup", [](const OpenMS::DeconvolvedSpectrum& self) -> const OpenMS::PeakGroup & { return self.getPrecursorPeakGroup(); }, nb::rv_policy::reference_internal, "Returns the precursor peak group (MSn, n>1)")
        .def("getPrecursorCharge", [](const OpenMS::DeconvolvedSpectrum& self) { return self.getPrecursorCharge(); }, "Returns the precursor charge")
        .def("getPrecursor", [](const OpenMS::DeconvolvedSpectrum& self) -> const OpenMS::Precursor & { return self.getPrecursor(); }, nb::rv_policy::reference_internal, "Returns the precursor peak")
        .def("getCurrentMaxMass", [](const OpenMS::DeconvolvedSpectrum& self, double max_mass) { return self.getCurrentMaxMass(max_mass); }, "max_mass"_a, "Returns the current max mass")
        .def("getCurrentMinMass", [](const OpenMS::DeconvolvedSpectrum& self, double min_mass) { return self.getCurrentMinMass(min_mass); }, "min_mass"_a, "Returns the current min mass")
        .def("getCurrentMaxAbsCharge", [](const OpenMS::DeconvolvedSpectrum& self, int max_abs_charge) { return self.getCurrentMaxAbsCharge(max_abs_charge); }, "max_abs_charge"_a, "Returns the current max charge")
        .def("getScanNumber", [](const OpenMS::DeconvolvedSpectrum& self) { return self.getScanNumber(); }, "Returns the scan number")
        .def("getPrecursorScanNumber", [](const OpenMS::DeconvolvedSpectrum& self) { return self.getPrecursorScanNumber(); }, "Returns the precursor scan number")
        .def("getQuantities", [](const OpenMS::DeconvolvedSpectrum& self) { return self.getQuantities(); }, "Returns isobaric quantities")
        .def("setQuantities", [](OpenMS::DeconvolvedSpectrum& self, const OpenMS::FLASHHelperClasses::IsobaricQuantities& quantities) { return self.setQuantities(quantities); }, "quantities"_a, "Sets isobaric quantities")
        .def("setPrecursor", [](OpenMS::DeconvolvedSpectrum& self, const OpenMS::Precursor& precursor) { return self.setPrecursor(precursor); }, "precursor"_a, "Sets the precursor")
        .def("setPrecursorScanNumber", [](OpenMS::DeconvolvedSpectrum& self, int scan_number) { return self.setPrecursorScanNumber(scan_number); }, "scan_number"_a, "Sets the precursor scan number")
        .def("setPrecursorPeakGroup", [](OpenMS::DeconvolvedSpectrum& self, const OpenMS::PeakGroup& pg) { return self.setPrecursorPeakGroup(pg); }, "pg"_a, "Sets the precursor peak group")
        .def("setOriginalSpectrum", [](OpenMS::DeconvolvedSpectrum& self, const OpenMS::MSSpectrum& spec) { return self.setOriginalSpectrum(spec); }, "spec"_a, "Sets the original spectrum")
        .def("setPeakGroups", [](OpenMS::DeconvolvedSpectrum& self, std::vector<OpenMS::PeakGroup>& x) { return self.setPeakGroups(x); }, "x"_a, "Sets peak groups")
        .def("begin", [](const OpenMS::DeconvolvedSpectrum& self) { return self.begin(); })
        .def("end", [](const OpenMS::DeconvolvedSpectrum& self) { return self.end(); })
        .def("push_back", [](OpenMS::DeconvolvedSpectrum& self, const OpenMS::PeakGroup& pg) { return self.push_back(pg); }, "pg"_a, "Adds a peak group")
        .def("pop_back", [](OpenMS::DeconvolvedSpectrum& self) { return self.pop_back(); }, "Removes the last peak group")
        .def("size", [](const OpenMS::DeconvolvedSpectrum& self) { return self.size(); }, "Returns number of peak groups")
        .def("clear", [](OpenMS::DeconvolvedSpectrum& self) { return self.clear(); }, "Clears all peak groups")
        .def("reserve", [](OpenMS::DeconvolvedSpectrum& self, size_t n) { return self.reserve(n); }, "n"_a, "Reserves space for n peak groups")
        .def("empty", [](const OpenMS::DeconvolvedSpectrum& self) { return self.empty(); }, "Returns true if no peak groups")
        .def("isDecoy", [](const OpenMS::DeconvolvedSpectrum& self) { return self.isDecoy(); }, "Returns true if this is a decoy spectrum")
        .def("sort", [](OpenMS::DeconvolvedSpectrum& self) { return self.sort(); }, "Sorts peak groups by monoisotopic mass")
        .def("sortByQscore", [](OpenMS::DeconvolvedSpectrum& self) { return self.sortByQscore(); }, "Sorts peak groups by Qscore")
        .def(nb::self < nb::self)
        .def(nb::self > nb::self)
        .def(nb::self == nb::self)
        .def("__iter__", [](OpenMS::DeconvolvedSpectrum& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::DeconvolvedSpectrum>(), "DeconvolvedSpectrum_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::DeconvolvedSpectrum& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::DeconvolvedSpectrum& self, size_t i) -> OpenMS::PeakGroup & { 
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::DeconvolvedSpectrum& self, size_t i, const OpenMS::PeakGroup& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a)
        ;

    // -----------------------------------------------------------------------
    // ExtractionCoordinates
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates>(m, "ExtractionCoordinates", "void extract_value_tophat # -> uses iterators")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates &>())
        .def("__copy__", [](const OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates& self) { return OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates(self); })
        .def("__deepcopy__", [](const OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates& self, nb::dict) { return OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates(self); }, "memo"_a)
        .def_rw("mz", &OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates::mz)
        .def_rw("mz_precursor", &OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates::mz_precursor)
        .def_rw("rt_start", &OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates::rt_start)
        .def_rw("rt_end", &OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates::rt_end)
        .def_rw("ion_mobility", &OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates::ion_mobility)
        .def_rw("id", &OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates::id)
        ;

    // -----------------------------------------------------------------------
    // FIAMSScheduler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FIAMSScheduler>(m, "FIAMSScheduler", "Scheduler for FIA-MS data processing")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FIAMSScheduler &>())
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, bool>())
        .def("__copy__", [](const OpenMS::FIAMSScheduler& self) { return OpenMS::FIAMSScheduler(self); })
        .def("__deepcopy__", [](const OpenMS::FIAMSScheduler& self, nb::dict) { return OpenMS::FIAMSScheduler(self); }, "memo"_a)
        .def("run", [](OpenMS::FIAMSScheduler& self) { return self.run(); }, "Run the FIA-MS data analysis for the batch defined in the @filename_")
        .def("getBaseDir", [](OpenMS::FIAMSScheduler& self) { return self.getBaseDir(); }, "Returns the base directory for the relevant paths from the csv file")
        ;

    // -----------------------------------------------------------------------
    // FLASHHelperClasses
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHHelperClasses>(m, "FLASHHelperClasses", 
        R"doc(
Wrapper struct for all the structs needed by FLASHDeconv.
Contains: PrecalculatedAveragine, MassFeature, IsobaricQuantities, LogMzPeak
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHHelperClasses &>())
        .def("__copy__", [](const OpenMS::FLASHHelperClasses& self) { return OpenMS::FLASHHelperClasses(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHHelperClasses& self, nb::dict) { return OpenMS::FLASHHelperClasses(self); }, "memo"_a)
        .def_static("getLogMz", [](double mz, bool positive) { return OpenMS::FLASHHelperClasses::getLogMz(mz, positive); }, "mz"_a, "positive"_a)
        .def_static("getChargeMass", [](bool positive_ioniziation_mode) { return OpenMS::FLASHHelperClasses::getChargeMass(positive_ioniziation_mode); }, "positive_ioniziation_mode"_a, "Calculate log mz from mz")
        ;

    // -----------------------------------------------------------------------
    // FeatureMapping
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureMapping>(m, "FeatureMapping", "OpenMS class FeatureMapping")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureMapping &>())
        .def("__copy__", [](const OpenMS::FeatureMapping& self) { return OpenMS::FeatureMapping(self); })
        .def("__deepcopy__", [](const OpenMS::FeatureMapping& self, nb::dict) { return OpenMS::FeatureMapping(self); }, "memo"_a)
        .def_static("assignMS2IndexToFeature", [](const OpenMS::MSExperiment& spectra, const OpenMS::FeatureMapping::FeatureMappingInfo& fm_info, const double& precursor_mz_tolerance, const double& precursor_rt_tolerance, bool ppm) { return OpenMS::FeatureMapping::assignMS2IndexToFeature(spectra, fm_info, precursor_mz_tolerance, precursor_rt_tolerance, ppm); }, "spectra"_a, "fm_info"_a, "precursor_mz_tolerance"_a, "precursor_rt_tolerance"_a, "ppm"_a)
        ;

    // -----------------------------------------------------------------------
    // FeatureMapping_FeatureMappingInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureMapping::FeatureMappingInfo>(m, "FeatureMapping_FeatureMappingInfo", "OpenMS class FeatureMapping_FeatureMappingInfo")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureMapping::FeatureMappingInfo &>())
        .def("__copy__", [](const OpenMS::FeatureMapping::FeatureMappingInfo& self) { return OpenMS::FeatureMapping::FeatureMappingInfo(self); })
        .def("__deepcopy__", [](const OpenMS::FeatureMapping::FeatureMappingInfo& self, nb::dict) { return OpenMS::FeatureMapping::FeatureMappingInfo(self); }, "memo"_a)
        .def_rw("feature_maps", &OpenMS::FeatureMapping::FeatureMappingInfo::feature_maps)
        .def_rw("kd_tree", &OpenMS::FeatureMapping::FeatureMappingInfo::kd_tree)
        ;

    // -----------------------------------------------------------------------
    // FeatureMapping_FeatureToMs2Indices
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureMapping::FeatureToMs2Indices>(m, "FeatureMapping_FeatureToMs2Indices", "OpenMS class FeatureMapping_FeatureToMs2Indices")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureMapping::FeatureToMs2Indices &>())
        .def("__copy__", [](const OpenMS::FeatureMapping::FeatureToMs2Indices& self) { return OpenMS::FeatureMapping::FeatureToMs2Indices(self); })
        .def("__deepcopy__", [](const OpenMS::FeatureMapping::FeatureToMs2Indices& self, nb::dict) { return OpenMS::FeatureMapping::FeatureToMs2Indices(self); }, "memo"_a)
        .def_ro("assignedMS2", &OpenMS::FeatureMapping::FeatureToMs2Indices::assignedMS2)
        .def_rw("unassignedMS2", &OpenMS::FeatureMapping::FeatureToMs2Indices::unassignedMS2)
        ;

    // -----------------------------------------------------------------------
    // HyperScore
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::HyperScore>(m, "HyperScore", "OpenMS class HyperScore")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::HyperScore &>())
        .def("__copy__", [](const OpenMS::HyperScore& self) { return OpenMS::HyperScore(self); })
        .def("__deepcopy__", [](const OpenMS::HyperScore& self, nb::dict) { return OpenMS::HyperScore(self); }, "memo"_a)
        .def_static("compute", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::MSSpectrum& theo_spectrum) { return OpenMS::HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum); }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "exp_spectrum"_a, "theo_spectrum"_a)
        .def_static("compute", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::DataArrays::IntegerDataArray& exp_charges, const OpenMS::MSSpectrum& theo_spectrum, const OpenMS::DataArrays::IntegerDataArray& theo_charges) { return OpenMS::HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, exp_charges, theo_spectrum, theo_charges); }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "exp_spectrum"_a, "exp_charges"_a, "theo_spectrum"_a, "theo_charges"_a)
        .def_static("compute", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::DataArrays::IntegerDataArray& exp_charges, const OpenMS::MSSpectrum& theo_spectrum, const OpenMS::DataArrays::IntegerDataArray& theo_charges, std::vector<double> intensity_sum) {
            auto score = OpenMS::HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, exp_charges, theo_spectrum, theo_charges, intensity_sum);
            return nb::make_tuple(score, intensity_sum);
        }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "exp_spectrum"_a, "exp_charges"_a, "theo_spectrum"_a, "theo_charges"_a, "intensity_sum"_a)
        ;

    // -----------------------------------------------------------------------
    // IDConflictResolverAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDConflictResolverAlgorithm>(m, "IDConflictResolverAlgorithm", "OpenMS class IDConflictResolverAlgorithm")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IDConflictResolverAlgorithm &>())
        .def("__copy__", [](const OpenMS::IDConflictResolverAlgorithm& self) { return OpenMS::IDConflictResolverAlgorithm(self); })
        .def("__deepcopy__", [](const OpenMS::IDConflictResolverAlgorithm& self, nb::dict) { return OpenMS::IDConflictResolverAlgorithm(self); }, "memo"_a)
        .def_static("resolve", [](OpenMS::FeatureMap& features, bool keep_matching) { return OpenMS::IDConflictResolverAlgorithm::resolve(features, keep_matching); }, "features"_a, "keep_matching"_a)
        .def_static("resolve", [](OpenMS::ConsensusMap& features, bool keep_matching) { return OpenMS::IDConflictResolverAlgorithm::resolve(features, keep_matching); }, "features"_a, "keep_matching"_a)
        .def_static("resolveAllHitRankAggregation", [](OpenMS::FeatureMap& features) { return OpenMS::IDConflictResolverAlgorithm::resolveAllHitRankAggregation(features); }, "features"_a,
            R"doc(
Resolves ambiguous annotations of features with peptide identifications using rank aggregation.
For each feature, peptide hits across all identifications are aggregated by rank.
Each unique sequence is assigned a rank in every identification in which it appears
(rank 0 = best hit). Sequences absent from an identification receive a penalty rank
equal to the maximum number of considered hits. The sequence with the best average
rank score is selected as the winner.
:param features: FeatureMap to work on (modified in-place)
)doc")
        .def_static("resolveAllHitRankAggregation", [](OpenMS::ConsensusMap& features) { return OpenMS::IDConflictResolverAlgorithm::resolveAllHitRankAggregation(features); }, "features"_a,
            R"doc(
Resolves ambiguous annotations of consensus features with peptide identifications using rank aggregation.
For each consensus feature, peptide hits across all identifications are aggregated by rank.
Each unique sequence is assigned a rank in every identification in which it appears
(rank 0 = best hit). Sequences absent from an identification receive a penalty rank
equal to the maximum number of considered hits. The sequence with the best average
rank score is selected as the winner.
:param features: ConsensusMap to work on (modified in-place)
)doc")
        .def_static("resolveBetweenFeatures", [](OpenMS::FeatureMap& features) { return OpenMS::IDConflictResolverAlgorithm::resolveBetweenFeatures(features); }, "features"_a, 
            R"doc(
Resolves ambiguous annotations of consensus features with peptide identifications\n
The the filtered identifications are added to the vector of unassigned peptides
and also reduced to a single best hit
:param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
)doc")
        .def_static("resolveBetweenFeatures", [](OpenMS::ConsensusMap& features) { return OpenMS::IDConflictResolverAlgorithm::resolveBetweenFeatures(features); }, "features"_a, 
            R"doc(
Resolves ambiguous annotations of consensus features with peptide identifications\n
The the filtered identifications are added to the vector of unassigned peptides
and also reduced to a single best hit
:param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
)doc")
        ;

    // -----------------------------------------------------------------------
    // IDMapper_PeptideIdentificationListState
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDMapper::PeptideIdentificationListState>(m, "IDMapper_PeptideIdentificationListState", ":return: A struct of vectors holding spectra indices of the partitioning")
        .def(nb::init<>())
        .def_rw("no_precursors", &OpenMS::IDMapper::PeptideIdentificationListState::no_precursors)
        .def_rw("identified", &OpenMS::IDMapper::PeptideIdentificationListState::identified)
        .def_rw("unidentified", &OpenMS::IDMapper::PeptideIdentificationListState::unidentified)
        ;

    // -----------------------------------------------------------------------
    // IonIdentityMolecularNetworking
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IonIdentityMolecularNetworking>(m, "IonIdentityMolecularNetworking", "Includes the necessary functions to generate filed required for GNPS ion identity molecular networking (IIMN).")
        .def(nb::init<>())
        .def_static("annotateConsensusMap", [](OpenMS::ConsensusMap& consensus_map) { return OpenMS::IonIdentityMolecularNetworking::annotateConsensusMap(consensus_map); }, "consensus_map"_a)
        .def_static("writeSupplementaryPairTable", [](const OpenMS::ConsensusMap& consensus_map, const OpenMS::String& output_file) { return OpenMS::IonIdentityMolecularNetworking::writeSupplementaryPairTable(consensus_map, output_file); }, "consensus_map"_a, "output_file"_a, 
            R"doc(
Annotate ConsensusMap for ion identity molecular networking (IIMN) workflow by GNPS.
Adds meta values Constants::UserParams::IIMN_ROW_ID (unique index for each feature), Constants::UserParams::IIMN_ADDUCT_PARTNERS (related features row IDs)
and Constants::UserParams::IIMN_ANNOTATION_NETWORK_NUMBER (all related features with different adduct states) get the same network number).
This method requires the features annotated with the Constants::UserParams::IIMN_LINKED_GROUPS meta value.
If at least one of the features has an annotation for Constants::UserParam::IIMN_LINKED_GROUPS, annotate ConsensusMap for IIMN.
:param consensus_map: Input ConsensusMap without IIMN annotations.
)doc")
        ;

    // -----------------------------------------------------------------------
    // IsobaricIsotopeCorrector
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsobaricIsotopeCorrector>(m, "IsobaricIsotopeCorrector", 
        R"doc(
Performs isotope impurity correction on intensities extracted from
isobaric labeling experiments
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IsobaricIsotopeCorrector &>())
        .def("__copy__", [](const OpenMS::IsobaricIsotopeCorrector& self) { return OpenMS::IsobaricIsotopeCorrector(self); })
        .def("__deepcopy__", [](const OpenMS::IsobaricIsotopeCorrector& self, nb::dict) { return OpenMS::IsobaricIsotopeCorrector(self); }, "memo"_a)
        .def_static("correctIsotopicImpurities", [](const OpenMS::ConsensusMap& consensus_map_in, const OpenMS::IsobaricQuantitationMethod* quant_method) {
            OpenMS::ConsensusMap consensus_map_out;
            auto stats = OpenMS::IsobaricIsotopeCorrector::correctIsotopicImpurities(consensus_map_in, consensus_map_out, quant_method);
            return std::make_pair(consensus_map_out, stats);
        }, "consensus_map_in"_a, "quant_method"_a, "Correct isotopic impurities in a ConsensusMap, returns (corrected_map, statistics)")
        ;

    // -----------------------------------------------------------------------
    // IsobaricQuantifierStatistics
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsobaricQuantifierStatistics>(m, "IsobaricQuantifierStatistics", 
        R"doc(
Statistics for isobaric quantitation performance. Tracks the number of
MS2 spectra processed, empty channels, and comparison metrics between
NNLS (non-negative least squares) and naive matrix inversion methods
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IsobaricQuantifierStatistics &>())
        .def("__copy__", [](const OpenMS::IsobaricQuantifierStatistics& self) { return OpenMS::IsobaricQuantifierStatistics(self); })
        .def("__deepcopy__", [](const OpenMS::IsobaricQuantifierStatistics& self, nb::dict) { return OpenMS::IsobaricQuantifierStatistics(self); }, "memo"_a)
        .def("reset", [](OpenMS::IsobaricQuantifierStatistics& self) { return self.reset(); }, "Resets all statistics counters to zero")
        .def_rw("channel_count", &OpenMS::IsobaricQuantifierStatistics::channel_count)
        .def_rw("iso_number_ms2_negative", &OpenMS::IsobaricQuantifierStatistics::iso_number_ms2_negative)
        .def_rw("iso_number_reporter_negative", &OpenMS::IsobaricQuantifierStatistics::iso_number_reporter_negative)
        .def_rw("iso_number_reporter_different", &OpenMS::IsobaricQuantifierStatistics::iso_number_reporter_different)
        .def_rw("iso_solution_different_intensity", &OpenMS::IsobaricQuantifierStatistics::iso_solution_different_intensity)
        .def_rw("iso_total_intensity_negative", &OpenMS::IsobaricQuantifierStatistics::iso_total_intensity_negative)
        .def_rw("number_ms2_total", &OpenMS::IsobaricQuantifierStatistics::number_ms2_total)
        .def_rw("number_ms2_empty", &OpenMS::IsobaricQuantifierStatistics::number_ms2_empty)
        .def_rw("empty_channels", &OpenMS::IsobaricQuantifierStatistics::empty_channels)
        ;

    // -----------------------------------------------------------------------
    // IsobaricQuantities
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHHelperClasses::IsobaricQuantities>(m, "IsobaricQuantities", 
        R"doc(
Isobaric quantities from isobaric quantification.
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHHelperClasses::IsobaricQuantities &>())
        .def("__copy__", [](const OpenMS::FLASHHelperClasses::IsobaricQuantities& self) { return OpenMS::FLASHHelperClasses::IsobaricQuantities(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHHelperClasses::IsobaricQuantities& self, nb::dict) { return OpenMS::FLASHHelperClasses::IsobaricQuantities(self); }, "memo"_a)
        .def("empty", [](const OpenMS::FLASHHelperClasses::IsobaricQuantities& self) { return self.empty(); }, "Returns true if no quantities stored")
        .def_rw("scan", &OpenMS::FLASHHelperClasses::IsobaricQuantities::scan)
        .def_rw("rt", &OpenMS::FLASHHelperClasses::IsobaricQuantities::rt)
        .def_rw("precursor_mz", &OpenMS::FLASHHelperClasses::IsobaricQuantities::precursor_mz)
        .def_rw("precursor_mass", &OpenMS::FLASHHelperClasses::IsobaricQuantities::precursor_mass)
        .def_rw("quantities", &OpenMS::FLASHHelperClasses::IsobaricQuantities::quantities)
        .def_rw("merged_quantities", &OpenMS::FLASHHelperClasses::IsobaricQuantities::merged_quantities)
        ;

    // -----------------------------------------------------------------------
    // ItraqConstants
    // -----------------------------------------------------------------------
    auto itraqconstants_class = nb::class_<OpenMS::ItraqConstants>(m, "ItraqConstants", 
        R"doc(
Some constants used throughout iTRAQ classes
Constants for iTRAQ experiments and a ChannelInfo structure to store information about a single channel
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ItraqConstants &>())
        .def("__copy__", [](const OpenMS::ItraqConstants& self) { return OpenMS::ItraqConstants(self); })
        .def("__deepcopy__", [](const OpenMS::ItraqConstants& self, nb::dict) { return OpenMS::ItraqConstants(self); }, "memo"_a)
        .def_static("getIsotopeMatrixAsStringList", [](int itraq_type, const std::vector<OpenMS::Matrix<double>>& isotope_corrections) { return OpenMS::ItraqConstants::getIsotopeMatrixAsStringList(itraq_type, isotope_corrections); }, "itraq_type"_a, "isotope_corrections"_a)
        .def_static("updateIsotopeMatrixFromStringList", [](int itraq_type, const std::vector<OpenMS::String>& channels, std::vector<OpenMS::Matrix<double>> isotope_corrections) {
            OpenMS::ItraqConstants::updateIsotopeMatrixFromStringList(itraq_type, channels, isotope_corrections);
            return isotope_corrections;
        }, "itraq_type"_a, "channels"_a, "isotope_corrections"_a,
            R"doc(
Convert isotope correction matrix to stringlist\n
Each line is converted into a string of the format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
Useful for creating parameters or debug output
:param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
:param isotope_corrections: Vector of the two matrices (4plex, 8plex)
:returns: Updated isotope_corrections
)doc")
        .def_static("translateIsotopeMatrix", [](const int& itraq_type, const std::vector<OpenMS::Matrix<double>>& isotope_corrections) { return OpenMS::ItraqConstants::translateIsotopeMatrix(itraq_type, isotope_corrections); }, "itraq_type"_a, "isotope_corrections"_a, 
            R"doc(
Convert strings to isotope correction matrix rows\n
Each string of format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
is parsed and the corresponding channel(row) in the matrix is updated
Not all channels need to be present, missing channels will be left untouched
Useful to update the matrix with user isotope correction values
:param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
:param channels: New channel isotope values as strings
:param isotope_corrections: Vector of the two matrices (4plex, 8plex)
)doc")
        ;
    // ITRAQ_TYPES enum nested under ItraqConstants
    auto itraq_types_enum = nb::enum_<OpenMS::ItraqConstants::ITRAQ_TYPES>(itraqconstants_class, "ITRAQ_TYPES", "Enum for iTRAQ/TMT channel types", nb::is_arithmetic())
        .value("FOURPLEX", OpenMS::ItraqConstants::ITRAQ_TYPES::FOURPLEX)
        .value("EIGHTPLEX", OpenMS::ItraqConstants::ITRAQ_TYPES::EIGHTPLEX)
        .value("TMT_SIXPLEX", OpenMS::ItraqConstants::ITRAQ_TYPES::TMT_SIXPLEX)
        .value("SIZE_OF_ITRAQ_TYPES", OpenMS::ItraqConstants::ITRAQ_TYPES::SIZE_OF_ITRAQ_TYPES)

        .export_values();
    // Module-level alias so pyopenms.ITRAQ_TYPES works
    m.attr("ITRAQ_TYPES") = itraq_types_enum;

    // -----------------------------------------------------------------------
    // ChannelInfo (ItraqConstants::ChannelInfo)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ItraqConstants::ChannelInfo>(m, "ChannelInfo",
        "Information about an iTRAQ/TMT channel")
        .def(nb::init<>())
        .def_rw("description", &OpenMS::ItraqConstants::ChannelInfo::description)
        .def_rw("name", &OpenMS::ItraqConstants::ChannelInfo::name)
        .def_rw("id", &OpenMS::ItraqConstants::ChannelInfo::id)
        .def_rw("center", &OpenMS::ItraqConstants::ChannelInfo::center)
        .def_rw("active", &OpenMS::ItraqConstants::ChannelInfo::active)
        ;

    // -----------------------------------------------------------------------
    // KDTreeFeatureNode
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::KDTreeFeatureNode>(m, "KDTreeFeatureNode", 
        R"doc(
A node of the kD-tree with pointer to corresponding data and index
)doc")
        .def(nb::init<const OpenMS::KDTreeFeatureNode &>())
        .def("__copy__", [](const OpenMS::KDTreeFeatureNode& self) { return OpenMS::KDTreeFeatureNode(self); })
        .def("__deepcopy__", [](const OpenMS::KDTreeFeatureNode& self, nb::dict) { return OpenMS::KDTreeFeatureNode(self); }, "memo"_a)
        .def("__getitem__", [](OpenMS::KDTreeFeatureNode& self, size_t i) { return self[i]; })
        .def("getIndex", [](const OpenMS::KDTreeFeatureNode& self) { return self.getIndex(); }, "Returns index of corresponding feature in data_")
        ;

    // -----------------------------------------------------------------------
    // LightCompound
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::LightCompound>(m, "LightCompound", "OpenMS class LightCompound")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::LightCompound &>())
        .def("__copy__", [](const OpenSwath::LightCompound& self) { return OpenSwath::LightCompound(self); })
        .def("__deepcopy__", [](const OpenSwath::LightCompound& self, nb::dict) { return OpenSwath::LightCompound(self); }, "memo"_a)
        .def("setDriftTime", [](OpenSwath::LightCompound& self, double d) { return self.setDriftTime(d); }, "d"_a)
        .def("getDriftTime", [](const OpenSwath::LightCompound& self) { return self.getDriftTime(); })
        .def("getChargeState", [](const OpenSwath::LightCompound& self) { return self.getChargeState(); })
        .def("isPeptide", [](const OpenSwath::LightCompound& self) { return self.isPeptide(); })
        .def("setChargeState", [](OpenSwath::LightCompound& self, int ch) { return self.setChargeState(ch); }, "ch"_a)
        .def_rw("rt", &OpenSwath::LightCompound::rt)
        .def_rw("drift_time", &OpenSwath::LightCompound::drift_time)
        .def_rw("charge", &OpenSwath::LightCompound::charge)
        .def_rw("sequence", &OpenSwath::LightCompound::sequence)
        .def_rw("protein_refs", &OpenSwath::LightCompound::protein_refs)
        .def_rw("peptide_group_label", &OpenSwath::LightCompound::peptide_group_label)
        .def_rw("gene_name", &OpenSwath::LightCompound::gene_name)
        .def_rw("id", &OpenSwath::LightCompound::id)
        .def_rw("sum_formula", &OpenSwath::LightCompound::sum_formula)
        .def_rw("compound_name", &OpenSwath::LightCompound::compound_name)
        .def_rw("label_type", &OpenSwath::LightCompound::label_type)
        .def_rw("smiles", &OpenSwath::LightCompound::smiles)
        .def_rw("adducts", &OpenSwath::LightCompound::adducts)
        .def_rw("modifications", &OpenSwath::LightCompound::modifications)
        ;

    // -----------------------------------------------------------------------
    // LightModification
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::LightModification>(m, "LightModification", "OpenMS class LightModification")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::LightModification &>())
        .def("__copy__", [](const OpenSwath::LightModification& self) { return OpenSwath::LightModification(self); })
        .def("__deepcopy__", [](const OpenSwath::LightModification& self, nb::dict) { return OpenSwath::LightModification(self); }, "memo"_a)
        .def_rw("location", &OpenSwath::LightModification::location)
        .def_rw("unimod_id", &OpenSwath::LightModification::unimod_id)
        ;

    // -----------------------------------------------------------------------
    // LightProtein
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::LightProtein>(m, "LightProtein", "OpenMS class LightProtein")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::LightProtein &>())
        .def("__copy__", [](const OpenSwath::LightProtein& self) { return OpenSwath::LightProtein(self); })
        .def("__deepcopy__", [](const OpenSwath::LightProtein& self, nb::dict) { return OpenSwath::LightProtein(self); }, "memo"_a)
        .def_rw("id", &OpenSwath::LightProtein::id)
        .def_rw("sequence", &OpenSwath::LightProtein::sequence)
        .def_rw("uniprot_id", &OpenSwath::LightProtein::uniprot_id)
        ;

    // -----------------------------------------------------------------------
    // LightTargetedExperiment
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::LightTargetedExperiment>(m, "LightTargetedExperiment", "OpenMS class LightTargetedExperiment")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::LightTargetedExperiment &>())
        .def("__copy__", [](const OpenSwath::LightTargetedExperiment& self) { return OpenSwath::LightTargetedExperiment(self); })
        .def("__deepcopy__", [](const OpenSwath::LightTargetedExperiment& self, nb::dict) { return OpenSwath::LightTargetedExperiment(self); }, "memo"_a)
        .def("getTransitions", [](const OpenSwath::LightTargetedExperiment& self) { return self.getTransitions(); })
        .def("getCompounds", [](const OpenSwath::LightTargetedExperiment& self) { return self.getCompounds(); })
        .def("getProteins", [](const OpenSwath::LightTargetedExperiment& self) { return self.getProteins(); })
        .def_rw("transitions", &OpenSwath::LightTargetedExperiment::transitions)
        .def_rw("compounds", &OpenSwath::LightTargetedExperiment::compounds)
        .def_rw("proteins", &OpenSwath::LightTargetedExperiment::proteins)

        .def("getCompoundByRef", [](OpenSwath::LightTargetedExperiment& self, const std::string& ref) { return self.getCompoundByRef(ref); }, "ref"_a)

        .def("getPeptideByRef", [](OpenSwath::LightTargetedExperiment& self, const std::string& ref) { return self.getPeptideByRef(ref); }, "ref"_a)
        ;

    // -----------------------------------------------------------------------
    // LightTransition
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::LightTransition>(m, "LightTransition", "OpenMS class LightTransition")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::LightTransition &>())
        .def("__copy__", [](const OpenSwath::LightTransition& self) { return OpenSwath::LightTransition(self); })
        .def("__deepcopy__", [](const OpenSwath::LightTransition& self, nb::dict) { return OpenSwath::LightTransition(self); }, "memo"_a)
        .def("getProductChargeState", [](const OpenSwath::LightTransition& self) { return self.getProductChargeState(); })
        .def("isProductChargeStateSet", [](const OpenSwath::LightTransition& self) { return self.isProductChargeStateSet(); })
        .def("getNativeID", [](const OpenSwath::LightTransition& self) { return self.getNativeID(); })
        .def("getPeptideRef", [](const OpenSwath::LightTransition& self) { return self.getPeptideRef(); })
        .def("getLibraryIntensity", [](const OpenSwath::LightTransition& self) { return self.getLibraryIntensity(); })
        .def("setLibraryIntensity", [](OpenSwath::LightTransition& self, double l) { return self.setLibraryIntensity(l); }, "l"_a)
        .def("getProductMZ", [](const OpenSwath::LightTransition& self) { return self.getProductMZ(); })
        .def("getPrecursorMZ", [](const OpenSwath::LightTransition& self) { return self.getPrecursorMZ(); })
        .def("isPrecursorImSet", [](const OpenSwath::LightTransition& self) { return self.isPrecursorImSet(); })
        .def("getPrecursorIM", [](const OpenSwath::LightTransition& self) { return self.getPrecursorIM(); })
        .def("getCompoundRef", [](const OpenSwath::LightTransition& self) { return self.getCompoundRef(); })
        .def("getDecoy", [](const OpenSwath::LightTransition& self) { return self.getDecoy(); })
        .def("setDecoy", [](OpenSwath::LightTransition& self, bool d) { return self.setDecoy(d); }, "d"_a)
        .def("getFragmentType", [](const OpenSwath::LightTransition& self) { return self.getFragmentType(); })
        .def("setFragmentType", [](OpenSwath::LightTransition& self, std::string s) { return self.setFragmentType(s); }, "s"_a)
        .def("getAnnotation", [](const OpenSwath::LightTransition& self) { return self.getAnnotation(); })
        .def("setDetectingTransition", [](OpenSwath::LightTransition& self, bool d) { return self.setDetectingTransition(d); }, "d"_a)
        .def("isDetectingTransition", [](const OpenSwath::LightTransition& self) { return self.isDetectingTransition(); })
        .def("setQuantifyingTransition", [](OpenSwath::LightTransition& self, bool q) { return self.setQuantifyingTransition(q); }, "q"_a)
        .def("isQuantifyingTransition", [](const OpenSwath::LightTransition& self) { return self.isQuantifyingTransition(); })
        .def("setIdentifyingTransition", [](OpenSwath::LightTransition& self, bool i) { return self.setIdentifyingTransition(i); }, "i"_a)
        .def("isIdentifyingTransition", [](const OpenSwath::LightTransition& self) { return self.isIdentifyingTransition(); })
        .def_rw("transition_name", &OpenSwath::LightTransition::transition_name)
        .def_rw("peptide_ref", &OpenSwath::LightTransition::peptide_ref)
        .def_rw("library_intensity", &OpenSwath::LightTransition::library_intensity)
        .def_rw("product_mz", &OpenSwath::LightTransition::product_mz)
        .def_rw("precursor_mz", &OpenSwath::LightTransition::precursor_mz)
        .def_rw("precursor_im", &OpenSwath::LightTransition::precursor_im)
        .def_rw("fragment_charge", &OpenSwath::LightTransition::fragment_charge)
        .def_rw("fragment_nr", &OpenSwath::LightTransition::fragment_nr)
        .def_rw("peptidoforms", &OpenSwath::LightTransition::peptidoforms)
        .def_prop_rw("decoy",
            [](const OpenSwath::LightTransition& t) { return t.getDecoy(); },
            [](OpenSwath::LightTransition& t, bool v) { t.setDecoy(v); })
        .def_prop_rw("detecting_transition",
            [](const OpenSwath::LightTransition& t) { return t.isDetectingTransition(); },
            [](OpenSwath::LightTransition& t, bool v) { t.setDetectingTransition(v); })
        .def_prop_rw("identifying_transition",
            [](const OpenSwath::LightTransition& t) { return t.isIdentifyingTransition(); },
            [](OpenSwath::LightTransition& t, bool v) { t.setIdentifyingTransition(v); })
        .def_prop_rw("quantifying_transition",
            [](const OpenSwath::LightTransition& t) { return t.isQuantifyingTransition(); },
            [](OpenSwath::LightTransition& t, bool v) { t.setQuantifyingTransition(v); })
        ;

    // -----------------------------------------------------------------------
    // LogMzPeak
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHHelperClasses::LogMzPeak>(m, "LogMzPeak", 
        R"doc(
Log transformed peak from original peak.
Contains information such as charge, isotope index, and uncharged mass.
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHHelperClasses::LogMzPeak &>())
        .def("__copy__", [](const OpenMS::FLASHHelperClasses::LogMzPeak& self) { return OpenMS::FLASHHelperClasses::LogMzPeak(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHHelperClasses::LogMzPeak& self, nb::dict) { return OpenMS::FLASHHelperClasses::LogMzPeak(self); }, "memo"_a)
        .def("getUnchargedMass", [](const OpenMS::FLASHHelperClasses::LogMzPeak& self) { return self.getUnchargedMass(); }, "Get uncharged mass (0 if no charge set)")
        .def(nb::self < nb::self)
        .def(nb::self > nb::self)
        .def(nb::self == nb::self)
        .def_rw("mz", &OpenMS::FLASHHelperClasses::LogMzPeak::mz)
        .def_rw("intensity", &OpenMS::FLASHHelperClasses::LogMzPeak::intensity)
        .def_rw("logMz", &OpenMS::FLASHHelperClasses::LogMzPeak::logMz)
        .def_rw("mass", &OpenMS::FLASHHelperClasses::LogMzPeak::mass)
        .def_rw("abs_charge", &OpenMS::FLASHHelperClasses::LogMzPeak::abs_charge)
        .def_rw("is_positive", &OpenMS::FLASHHelperClasses::LogMzPeak::is_positive)
        .def_rw("isotopeIndex", &OpenMS::FLASHHelperClasses::LogMzPeak::isotopeIndex)

        .def("__init__", [](OpenMS::FLASHHelperClasses::LogMzPeak* self, const OpenMS::Peak1D& peak, bool positive) {
            new (self) OpenMS::FLASHHelperClasses::LogMzPeak(peak, positive);
        }, "peak"_a, "positive"_a)
        ;

    // -----------------------------------------------------------------------
    // MRMBatchFeatureSelector
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMBatchFeatureSelector>(m, "MRMBatchFeatureSelector", "OpenMS class MRMBatchFeatureSelector")
        ;

    // -----------------------------------------------------------------------
    // MRMFQC_ComponentGroupPairQCs
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureQC::ComponentGroupPairQCs>(m, "MRMFQC_ComponentGroupPairQCs", "OpenMS class MRMFQC_ComponentGroupPairQCs")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeatureQC::ComponentGroupPairQCs &>())
        .def("__copy__", [](const OpenMS::MRMFeatureQC::ComponentGroupPairQCs& self) { return OpenMS::MRMFeatureQC::ComponentGroupPairQCs(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeatureQC::ComponentGroupPairQCs& self, nb::dict) { return OpenMS::MRMFeatureQC::ComponentGroupPairQCs(self); }, "memo"_a)
        .def_rw("component_group_name", &OpenMS::MRMFeatureQC::ComponentGroupPairQCs::component_group_name)
        .def_rw("resolution_pair_name", &OpenMS::MRMFeatureQC::ComponentGroupPairQCs::resolution_pair_name)
        .def_rw("resolution_l", &OpenMS::MRMFeatureQC::ComponentGroupPairQCs::resolution_l)
        .def_rw("resolution_u", &OpenMS::MRMFeatureQC::ComponentGroupPairQCs::resolution_u)
        .def_rw("rt_diff_l", &OpenMS::MRMFeatureQC::ComponentGroupPairQCs::rt_diff_l)
        .def_rw("rt_diff_u", &OpenMS::MRMFeatureQC::ComponentGroupPairQCs::rt_diff_u)
        ;

    // -----------------------------------------------------------------------
    // MRMFQC_ComponentGroupQCs
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureQC::ComponentGroupQCs>(m, "MRMFQC_ComponentGroupQCs", "OpenMS class MRMFQC_ComponentGroupQCs")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeatureQC::ComponentGroupQCs &>())
        .def("__copy__", [](const OpenMS::MRMFeatureQC::ComponentGroupQCs& self) { return OpenMS::MRMFeatureQC::ComponentGroupQCs(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeatureQC::ComponentGroupQCs& self, nb::dict) { return OpenMS::MRMFeatureQC::ComponentGroupQCs(self); }, "memo"_a)
        .def_rw("component_group_name", &OpenMS::MRMFeatureQC::ComponentGroupQCs::component_group_name)
        .def_rw("retention_time_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::retention_time_l)
        .def_rw("retention_time_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::retention_time_u)
        .def_rw("intensity_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::intensity_l)
        .def_rw("intensity_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::intensity_u)
        .def_rw("overall_quality_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::overall_quality_l)
        .def_rw("overall_quality_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::overall_quality_u)
        .def_rw("n_heavy_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_heavy_l)
        .def_rw("n_heavy_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_heavy_u)
        .def_rw("n_light_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_light_l)
        .def_rw("n_light_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_light_u)
        .def_rw("n_detecting_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_detecting_l)
        .def_rw("n_detecting_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_detecting_u)
        .def_rw("n_quantifying_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_quantifying_l)
        .def_rw("n_quantifying_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_quantifying_u)
        .def_rw("n_identifying_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_identifying_l)
        .def_rw("n_identifying_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_identifying_u)
        .def_rw("n_transitions_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_transitions_l)
        .def_rw("n_transitions_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::n_transitions_u)
        .def_rw("ion_ratio_pair_name_1", &OpenMS::MRMFeatureQC::ComponentGroupQCs::ion_ratio_pair_name_1)
        .def_rw("ion_ratio_pair_name_2", &OpenMS::MRMFeatureQC::ComponentGroupQCs::ion_ratio_pair_name_2)
        .def_rw("ion_ratio_l", &OpenMS::MRMFeatureQC::ComponentGroupQCs::ion_ratio_l)
        .def_rw("ion_ratio_u", &OpenMS::MRMFeatureQC::ComponentGroupQCs::ion_ratio_u)
        .def_rw("ion_ratio_feature_name", &OpenMS::MRMFeatureQC::ComponentGroupQCs::ion_ratio_feature_name)
        ;

    // -----------------------------------------------------------------------
    // MRMFQC_ComponentQCs
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureQC::ComponentQCs>(m, "MRMFQC_ComponentQCs", "OpenMS class MRMFQC_ComponentQCs")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeatureQC::ComponentQCs &>())
        .def("__copy__", [](const OpenMS::MRMFeatureQC::ComponentQCs& self) { return OpenMS::MRMFeatureQC::ComponentQCs(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeatureQC::ComponentQCs& self, nb::dict) { return OpenMS::MRMFeatureQC::ComponentQCs(self); }, "memo"_a)
        .def_rw("component_name", &OpenMS::MRMFeatureQC::ComponentQCs::component_name)
        .def_rw("retention_time_l", &OpenMS::MRMFeatureQC::ComponentQCs::retention_time_l)
        .def_rw("retention_time_u", &OpenMS::MRMFeatureQC::ComponentQCs::retention_time_u)
        .def_rw("intensity_l", &OpenMS::MRMFeatureQC::ComponentQCs::intensity_l)
        .def_rw("intensity_u", &OpenMS::MRMFeatureQC::ComponentQCs::intensity_u)
        .def_rw("overall_quality_l", &OpenMS::MRMFeatureQC::ComponentQCs::overall_quality_l)
        .def_rw("overall_quality_u", &OpenMS::MRMFeatureQC::ComponentQCs::overall_quality_u)
        ;

    // -----------------------------------------------------------------------
    // MRMFeaturePicker
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeaturePicker>(m, "MRMFeaturePicker",
        R"doc(
_MRMFeaturePicker_ defines the structures containing parameters to be
used in [MRMTransitionGroupPicker](@ref MRMTransitionGroupPicker) for
components and components groups
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MRMFeaturePicker& self) { return OpenMS::MRMFeaturePicker(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeaturePicker& self, nb::dict) { return OpenMS::MRMFeaturePicker(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // MRMFP_ComponentParams (MRMFeaturePicker::ComponentParams)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeaturePicker::ComponentParams>(m, "MRMFP_ComponentParams",
        "Parameters for a single MRM component")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MRMFeaturePicker::ComponentParams& self) { return OpenMS::MRMFeaturePicker::ComponentParams(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeaturePicker::ComponentParams& self, nb::dict) { return OpenMS::MRMFeaturePicker::ComponentParams(self); }, "memo"_a)
        .def_rw("component_name", &OpenMS::MRMFeaturePicker::ComponentParams::component_name)
        .def_rw("component_group_name", &OpenMS::MRMFeaturePicker::ComponentParams::component_group_name)
        .def_rw("params", &OpenMS::MRMFeaturePicker::ComponentParams::params)
        ;

    // -----------------------------------------------------------------------
    // MRMFP_ComponentGroupParams (MRMFeaturePicker::ComponentGroupParams)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeaturePicker::ComponentGroupParams>(m, "MRMFP_ComponentGroupParams",
        "Parameters for a component group in MRM feature picking")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MRMFeaturePicker::ComponentGroupParams& self) { return OpenMS::MRMFeaturePicker::ComponentGroupParams(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeaturePicker::ComponentGroupParams& self, nb::dict) { return OpenMS::MRMFeaturePicker::ComponentGroupParams(self); }, "memo"_a)
        .def_rw("component_group_name", &OpenMS::MRMFeaturePicker::ComponentGroupParams::component_group_name)
        .def_rw("params", &OpenMS::MRMFeaturePicker::ComponentGroupParams::params)
        ;

    // -----------------------------------------------------------------------
    // MRMFeatureQC
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureQC>(m, "MRMFeatureQC",
        R"doc(
The MRMFeatureQC is a class to handle the parameters and options for
MRMFeatureFilter
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MRMFeatureQC& self) { return OpenMS::MRMFeatureQC(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeatureQC& self, nb::dict) { return OpenMS::MRMFeatureQC(self); }, "memo"_a)

        .def_prop_rw("component_qcs",
            [](OpenMS::MRMFeatureQC& self) -> std::vector<OpenMS::MRMFeatureQC::ComponentQCs>& { return self.component_qcs; },
            [](OpenMS::MRMFeatureQC& self, std::vector<OpenMS::MRMFeatureQC::ComponentQCs> v) { self.component_qcs = std::move(v); })

        .def_prop_rw("component_group_qcs",
            [](OpenMS::MRMFeatureQC& self) -> std::vector<OpenMS::MRMFeatureQC::ComponentGroupQCs>& { return self.component_group_qcs; },
            [](OpenMS::MRMFeatureQC& self, std::vector<OpenMS::MRMFeatureQC::ComponentGroupQCs> v) { self.component_group_qcs = std::move(v); })

        .def_prop_rw("component_group_pair_qcs",
            [](OpenMS::MRMFeatureQC& self) -> std::vector<OpenMS::MRMFeatureQC::ComponentGroupPairQCs>& { return self.component_group_pair_qcs; },
            [](OpenMS::MRMFeatureQC& self, std::vector<OpenMS::MRMFeatureQC::ComponentGroupPairQCs> v) { self.component_group_pair_qcs = std::move(v); })
        ;

    // -----------------------------------------------------------------------
    // MRMFeatureSelectorQMIP
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureSelectorQMIP>(m, "MRMFeatureSelectorQMIP", 
        R"doc(
Select MRM features using Quadratic Mixed Integer Programming (QMIP)
This class selects MRMFeatures based on relative retention time using a
quadratic mixed integer programming formulation. The method optimizes
feature selection while maintaining retention time consistency between
neighboring transitions.
Example usage:
.. code-block:: python
from pyopenms import *
# Load features
features = FeatureMap()
FeatureXMLFile().load("features.featureXML", features)
# Configure parameters
params = SelectorParameters()
params.nn_threshold = 4
params.optimal_threshold = 0.5
# Run QMIP selection
selector = MRMFeatureSelectorQMIP()
selected = FeatureMap()
selector.selectMRMFeature(features, selected, params)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeatureSelectorQMIP &>())
        .def("__copy__", [](const OpenMS::MRMFeatureSelectorQMIP& self) { return OpenMS::MRMFeatureSelectorQMIP(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeatureSelectorQMIP& self, nb::dict) { return OpenMS::MRMFeatureSelectorQMIP(self); }, "memo"_a)
        .def("selectMRMFeature", [](const OpenMS::MRMFeatureSelectorQMIP& self, const OpenMS::FeatureMap& features, OpenMS::FeatureMap& selected_filtered, const OpenMS::MRMFeatureSelector::SelectorParameters& parameters) { return self.selectMRMFeature(features, selected_filtered, parameters); }, "features"_a, "selected_filtered"_a, "parameters"_a)
        ;

    // -----------------------------------------------------------------------
    // MRMFeatureSelectorScore
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureSelectorScore>(m, "MRMFeatureSelectorScore", 
        R"doc(
Select MRM features using score-weighted linear programming
This class selects MRMFeatures based on a linear programming formulation
where each possible transition is weighted by a user-defined score
(typically retention time and peak intensity). Multiple scoring functions
can be combined using the score_weights parameter.
Example usage:
.. code-block:: python
from pyopenms import *
# Load features
features = FeatureMap()
FeatureXMLFile().load("features.featureXML", features)
# Configure parameters
params = SelectorParameters()
params.nn_threshold = 4
params.optimal_threshold = 0.5
# Run score-based selection
selector = MRMFeatureSelectorScore()
selected = FeatureMap()
selector.selectMRMFeature(features, selected, params)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeatureSelectorScore &>())
        .def("__copy__", [](const OpenMS::MRMFeatureSelectorScore& self) { return OpenMS::MRMFeatureSelectorScore(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeatureSelectorScore& self, nb::dict) { return OpenMS::MRMFeatureSelectorScore(self); }, "memo"_a)
        .def("selectMRMFeature", [](const OpenMS::MRMFeatureSelectorScore& self, const OpenMS::FeatureMap& features, OpenMS::FeatureMap& selected_filtered, const OpenMS::MRMFeatureSelector::SelectorParameters& parameters) { return self.selectMRMFeature(features, selected_filtered, parameters); }, "features"_a, "selected_filtered"_a, "parameters"_a)
        ;

    // -----------------------------------------------------------------------
    // MRMIonSeries
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMIonSeries>(m, "MRMIonSeries",
        R"doc(
Generate theoretical fragment ion series for use in MRMAssay and
MRMDecoy
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MRMIonSeries& self) { return OpenMS::MRMIonSeries(self); })
        .def("__deepcopy__", [](const OpenMS::MRMIonSeries& self, nb::dict) { return OpenMS::MRMIonSeries(self); }, "memo"_a)
        .def("annotateTransitionCV", [](OpenMS::MRMIonSeries& self, OpenMS::ReactionMonitoringTransition& tr, const OpenMS::String& annotation) { return self.annotateTransitionCV(tr, annotation); }, "tr"_a, "annotation"_a)
        .def("annotateTransition", [](OpenMS::MRMIonSeries& self, OpenMS::ReactionMonitoringTransition& tr, const OpenMS::TargetedExperiment::Peptide& peptide, double precursor_mz_threshold, double product_mz_threshold, bool enable_reannotation, const std::vector<OpenMS::String>& fragment_types, const std::vector<size_t>& fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow) { self.annotateTransition(tr, peptide, precursor_mz_threshold, product_mz_threshold, enable_reannotation, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, round_decPow); }, "tr"_a, "peptide"_a, "precursor_mz_threshold"_a, "product_mz_threshold"_a, "enable_reannotation"_a, "fragment_types"_a, "fragment_charges"_a, "enable_specific_losses"_a, "enable_unspecific_losses"_a, "round_decPow"_a = -4, "Annotates a transition with the corresponding fragment ion")
        ;

    // -----------------------------------------------------------------------
    // MRMRTNormalizer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMRTNormalizer>(m, "MRMRTNormalizer", "The MRMRTNormalizer will find retention time peptides in data")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMRTNormalizer &>())
        .def("__copy__", [](const OpenMS::MRMRTNormalizer& self) { return OpenMS::MRMRTNormalizer(self); })
        .def("__deepcopy__", [](const OpenMS::MRMRTNormalizer& self, nb::dict) { return OpenMS::MRMRTNormalizer(self); }, "memo"_a)
        .def_static("removeOutliersRANSAC", [](const std::vector<std::pair<double, double>>& pairs, double rsq_limit, double coverage_limit, size_t max_iterations, double max_rt_threshold, size_t sampling_size) { return OpenMS::MRMRTNormalizer::removeOutliersRANSAC(pairs, rsq_limit, coverage_limit, max_iterations, max_rt_threshold, sampling_size); }, "pairs"_a, "rsq_limit"_a, "coverage_limit"_a, "max_iterations"_a, "max_rt_threshold"_a, "sampling_size"_a)
        .def_static("removeOutliersIterative", [](const std::vector<std::pair<double, double>>& pairs, double rsq_limit, double coverage_limit, bool use_chauvenet, const OpenMS::String& method) { return OpenMS::MRMRTNormalizer::removeOutliersIterative(pairs, rsq_limit, coverage_limit, use_chauvenet, method); }, "pairs"_a, "rsq_limit"_a, "coverage_limit"_a, "use_chauvenet"_a, "method"_a)
        .def_static("chauvenet_probability", [](const std::vector<double>& residuals, int pos) { return OpenMS::MRMRTNormalizer::chauvenet_probability(residuals, pos); }, "residuals"_a, "pos"_a)
        .def_static("chauvenet", [](const std::vector<double>& residuals, int pos) { return OpenMS::MRMRTNormalizer::chauvenet(residuals, pos); }, "residuals"_a, "pos"_a)
        .def_static("computeBinnedCoverage", [](const std::pair<double, double>& rtRange, const std::vector<std::pair<double, double>>& pairs, int nrBins, int minPeptidesPerBin, int minBinsFilled) { return OpenMS::MRMRTNormalizer::computeBinnedCoverage(rtRange, pairs, nrBins, minPeptidesPerBin, minBinsFilled); }, "rtRange"_a, "pairs"_a, "nrBins"_a, "minPeptidesPerBin"_a, "minBinsFilled"_a)
        ;

    // -----------------------------------------------------------------------
    // MRMScoring
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::MRMScoring>(m, "MRMScoring", "This class implements different scores for peaks found in SRM/MRM")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::MRMScoring &>())
        .def("__copy__", [](const OpenSwath::MRMScoring& self) { return OpenSwath::MRMScoring(self); })
        .def("__deepcopy__", [](const OpenSwath::MRMScoring& self, nb::dict) { return OpenSwath::MRMScoring(self); }, "memo"_a)
        .def("calcXcorrCoelutionScore", [](OpenSwath::MRMScoring& self) { return self.calcXcorrCoelutionScore(); }, "Calculate the cross-correlation coelution score. The score is a distance where zero indicates perfect coelution")
        .def("calcXcorrCoelutionWeightedScore", [](OpenSwath::MRMScoring& self, const std::vector<double>& normalized_library_intensity) { return self.calcXcorrCoelutionWeightedScore(normalized_library_intensity); }, "normalized_library_intensity"_a)
        .def("calcSeparateXcorrContrastCoelutionScore", [](OpenSwath::MRMScoring& self) { return self.calcSeparateXcorrContrastCoelutionScore(); }, "Calculate the separate cross-correlation contrast score")
        .def("calcXcorrPrecursorContrastCoelutionScore", [](OpenSwath::MRMScoring& self) { return self.calcXcorrPrecursorContrastCoelutionScore(); })
        .def("calcXcorrShapeScore", [](OpenSwath::MRMScoring& self) { return self.calcXcorrShapeScore(); }, 
            R"doc(
Calculate the precursor cross-correlation contrast score against the transitions
The score is a distance where zero indicates perfect coelution
)doc")
        .def("calcXcorrShapeWeightedScore", [](OpenSwath::MRMScoring& self, const std::vector<double>& normalized_library_intensity) { return self.calcXcorrShapeWeightedScore(normalized_library_intensity); }, "normalized_library_intensity"_a, 
            R"doc(
Calculate the cross-correlation shape score
The score is a correlation measure where 1 indicates perfect correlation
and 0 means no correlation.
)doc")
        .def("calcSeparateXcorrContrastShapeScore", [](OpenSwath::MRMScoring& self) { return self.calcSeparateXcorrContrastShapeScore(); }, "Calculate the separate cross-correlation contrast shape score")
        .def("calcXcorrPrecursorContrastShapeScore", [](OpenSwath::MRMScoring& self) { return self.calcXcorrPrecursorContrastShapeScore(); }, "Calculate the precursor cross-correlation shape score against the transitions")
        .def_static("calcRTScore", [](const OpenSwath::LightCompound& peptide, double normalized_experimental_rt) { return OpenSwath::MRMScoring::calcRTScore(peptide, normalized_experimental_rt); }, "peptide"_a, "normalized_experimental_rt"_a)
        .def("calcMIScore", [](OpenSwath::MRMScoring& self) { return self.calcMIScore(); })
        .def("calcMIWeightedScore", [](OpenSwath::MRMScoring& self, const std::vector<double>& normalized_library_intensity) { return self.calcMIWeightedScore(normalized_library_intensity); }, "normalized_library_intensity"_a)
        .def("calcMIPrecursorScore", [](OpenSwath::MRMScoring& self) { return self.calcMIPrecursorScore(); })
        .def("calcMIPrecursorContrastScore", [](OpenSwath::MRMScoring& self) { return self.calcMIPrecursorContrastScore(); })
        .def("calcMIPrecursorCombinedScore", [](OpenSwath::MRMScoring& self) { return self.calcMIPrecursorCombinedScore(); })
        .def("calcSeparateMIContrastScore", [](OpenSwath::MRMScoring& self) { return self.calcSeparateMIContrastScore(); })
        .def("getMIMatrix", [](const OpenSwath::MRMScoring& self) -> const OpenMS::Matrix<double>& { return self.getMIMatrix(); }, nb::rv_policy::reference_internal, "Returns the MI matrix")
        ;

    // -----------------------------------------------------------------------
    // MapAlignmentAlgorithmKD
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MapAlignmentAlgorithmKD>(m, "MapAlignmentAlgorithmKD", 
        R"doc(
An efficient reference-free feature map alignment algorithm for unlabeled data
This algorithm uses a kd-tree to efficiently compute conflict-free connected components (CCC)
in a compatibility graph on feature data. This graph is comprised of nodes corresponding
to features and edges connecting features f and f' iff both are within each other's tolerance
windows (wrt. RT and m/z difference). CCCs are those CCs that do not contain multiple features
from the same input map, and whose features all have the same charge state
All CCCs above a user-specified minimum size are considered true sets of corresponding features
and based on these, LOWESS transformations are computed for each input map such that the average
deviation from the mean retention time within all CCCs is minimized
private
)doc")
        .def(nb::init<size_t, OpenMS::Param>())
        .def("__copy__", [](const OpenMS::MapAlignmentAlgorithmKD& self) { return OpenMS::MapAlignmentAlgorithmKD(self); })
        .def("__deepcopy__", [](const OpenMS::MapAlignmentAlgorithmKD& self, nb::dict) { return OpenMS::MapAlignmentAlgorithmKD(self); }, "memo"_a)
        .def("addRTFitData", [](OpenMS::MapAlignmentAlgorithmKD& self, const OpenMS::KDTreeFeatureMaps& kd_data) { return self.addRTFitData(kd_data); }, "kd_data"_a, "Compute data points needed for RT transformation in the current `kd_data`, add to `fit_data_`")
        .def("fitLOWESS", [](OpenMS::MapAlignmentAlgorithmKD& self) { return self.fitLOWESS(); }, "Fit LOWESS to fit_data_, store final models in `transformations_`")
        .def("transform", [](const OpenMS::MapAlignmentAlgorithmKD& self, OpenMS::KDTreeFeatureMaps& kd_data) { return self.transform(kd_data); }, "kd_data"_a, "Transform RTs for `kd_data`")
        ;

    // -----------------------------------------------------------------------
    // MapAlignmentEvaluationAlgorithm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MapAlignmentEvaluationAlgorithm>(m, "MapAlignmentEvaluationAlgorithm", "Base class for all Caap evaluation algorithms")
        ;

    // -----------------------------------------------------------------------
    // MapAlignmentEvaluationAlgorithmPrecision
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MapAlignmentEvaluationAlgorithmPrecision, OpenMS::MapAlignmentEvaluationAlgorithm>(m, "MapAlignmentEvaluationAlgorithmPrecision", 
        R"doc(
Caap evaluation algorithm to obtain a precision value
MapAlignmentEvaluationAlgorithm
)doc")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // MapAlignmentEvaluationAlgorithmRecall
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MapAlignmentEvaluationAlgorithmRecall, OpenMS::MapAlignmentEvaluationAlgorithm>(m, "MapAlignmentEvaluationAlgorithmRecall", 
        R"doc(
Caap evaluation algorithm to obtain a recall value
MapAlignmentEvaluationAlgorithm
)doc")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // MapAlignmentTransformer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MapAlignmentTransformer>(m, "MapAlignmentTransformer", "This class collects functions for applying retention time transformations to data structures")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MapAlignmentTransformer &>())
        .def("__copy__", [](const OpenMS::MapAlignmentTransformer& self) { return OpenMS::MapAlignmentTransformer(self); })
        .def("__deepcopy__", [](const OpenMS::MapAlignmentTransformer& self, nb::dict) { return OpenMS::MapAlignmentTransformer(self); }, "memo"_a)
        .def_static("transformRetentionTimes", [](OpenMS::MSExperiment& msexp, const OpenMS::TransformationDescription& trafo, bool store_original_rt) { return OpenMS::MapAlignmentTransformer::transformRetentionTimes(msexp, trafo, store_original_rt); }, "msexp"_a, "trafo"_a, "store_original_rt"_a, "Applies the given transformation to a peak map")
        .def_static("transformRetentionTimes", [](OpenMS::FeatureMap& fmap, const OpenMS::TransformationDescription& trafo, bool store_original_rt) { return OpenMS::MapAlignmentTransformer::transformRetentionTimes(fmap, trafo, store_original_rt); }, "fmap"_a, "trafo"_a, "store_original_rt"_a, "Applies the given transformation to a peak map")
        .def_static("transformRetentionTimes", [](OpenMS::ConsensusMap& cmap, const OpenMS::TransformationDescription& trafo, bool store_original_rt) { return OpenMS::MapAlignmentTransformer::transformRetentionTimes(cmap, trafo, store_original_rt); }, "cmap"_a, "trafo"_a, "store_original_rt"_a, "Applies the given transformation to a peak map")
        .def_static("transformRetentionTimes", [](OpenMS::PeptideIdentificationList& pep_ids, const OpenMS::TransformationDescription& trafo, bool store_original_rt) { return OpenMS::MapAlignmentTransformer::transformRetentionTimes(pep_ids, trafo, store_original_rt); }, "pep_ids"_a, "trafo"_a, "store_original_rt"_a, "Applies the given transformation to a peak map")
        .def_static("transformRetentionTimes", [](OpenMS::IdentificationData& id_data, const OpenMS::TransformationDescription& trafo, bool store_original_rt) { return OpenMS::MapAlignmentTransformer::transformRetentionTimes(id_data, trafo, store_original_rt); }, "id_data"_a, "trafo"_a, "store_original_rt"_a, "Applies the given transformation to a peak map")
        ;

    // -----------------------------------------------------------------------
    // MassFeature_FDHS
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHHelperClasses::MassFeature>(m, "MassFeature_FDHS", 
        R"doc(
Mass feature (Deconvolved masses in spectra are traced to generate mass features).
Similar to LC-MS features but for deconvolved masses.
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHHelperClasses::MassFeature &>())
        .def("__copy__", [](const OpenMS::FLASHHelperClasses::MassFeature& self) { return OpenMS::FLASHHelperClasses::MassFeature(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHHelperClasses::MassFeature& self, nb::dict) { return OpenMS::FLASHHelperClasses::MassFeature(self); }, "memo"_a)
        .def(nb::self < nb::self)
        .def(nb::self > nb::self)
        .def(nb::self == nb::self)
        .def_rw("index", &OpenMS::FLASHHelperClasses::MassFeature::index)
        .def_rw("mt", &OpenMS::FLASHHelperClasses::MassFeature::mt)
        .def_rw("per_charge_intensity", &OpenMS::FLASHHelperClasses::MassFeature::per_charge_intensity)
        .def_rw("per_isotope_intensity", &OpenMS::FLASHHelperClasses::MassFeature::per_isotope_intensity)
        .def_rw("iso_offset", &OpenMS::FLASHHelperClasses::MassFeature::iso_offset)
        .def_rw("scan_number", &OpenMS::FLASHHelperClasses::MassFeature::scan_number)
        .def_rw("min_scan_number", &OpenMS::FLASHHelperClasses::MassFeature::min_scan_number)
        .def_rw("max_scan_number", &OpenMS::FLASHHelperClasses::MassFeature::max_scan_number)
        .def_rw("rep_charge", &OpenMS::FLASHHelperClasses::MassFeature::rep_charge)
        .def_rw("avg_mass", &OpenMS::FLASHHelperClasses::MassFeature::avg_mass)
        .def_rw("min_charge", &OpenMS::FLASHHelperClasses::MassFeature::min_charge)
        .def_rw("max_charge", &OpenMS::FLASHHelperClasses::MassFeature::max_charge)
        .def_rw("charge_count", &OpenMS::FLASHHelperClasses::MassFeature::charge_count)
        .def_rw("isotope_score", &OpenMS::FLASHHelperClasses::MassFeature::isotope_score)
        .def_rw("qscore", &OpenMS::FLASHHelperClasses::MassFeature::qscore)
        .def_rw("rep_mz", &OpenMS::FLASHHelperClasses::MassFeature::rep_mz)
        .def_rw("is_decoy", &OpenMS::FLASHHelperClasses::MassFeature::is_decoy)
        .def_rw("ms_level", &OpenMS::FLASHHelperClasses::MassFeature::ms_level)
        ;

    // -----------------------------------------------------------------------
    // MetaboTargetedAssay
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaboTargetedAssay>(m, "MetaboTargetedAssay", 
        R"doc(
This class provides methods for the extraction of targeted assays for
metabolomics
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaboTargetedAssay &>())
        .def("__copy__", [](const OpenMS::MetaboTargetedAssay& self) { return OpenMS::MetaboTargetedAssay(self); })
        .def("__deepcopy__", [](const OpenMS::MetaboTargetedAssay& self, nb::dict) { return OpenMS::MetaboTargetedAssay(self); }, "memo"_a)
        .def_rw("precursor_int", &OpenMS::MetaboTargetedAssay::precursor_int)
        .def_rw("transition_quality_score", &OpenMS::MetaboTargetedAssay::transition_quality_score)
        .def_rw("precursor_mz", &OpenMS::MetaboTargetedAssay::precursor_mz)
        .def_rw("compound_rt", &OpenMS::MetaboTargetedAssay::compound_rt)
        .def_rw("molecular_formula", &OpenMS::MetaboTargetedAssay::molecular_formula)
        .def_rw("compound_file", &OpenMS::MetaboTargetedAssay::compound_file)
        .def_rw("compound_name", &OpenMS::MetaboTargetedAssay::compound_name)
        .def_rw("compound_adduct", &OpenMS::MetaboTargetedAssay::compound_adduct)
        .def_rw("potential_cmp", &OpenMS::MetaboTargetedAssay::potential_cmp)
        .def_rw("potential_rmts", &OpenMS::MetaboTargetedAssay::potential_rmts)
        .def_static("extractMetaboTargetedAssay", [](const OpenMS::MSExperiment& spectra, const OpenMS::FeatureMapping::FeatureToMs2Indices& feature_ms2_index, const double& precursor_rt_tol, const double& precursor_mz_distance, const double& cosine_sim_threshold, const double& transition_threshold, const double& min_fragment_mz, const double& max_fragment_mz, const bool& method_consensus_spectrum, const bool& exclude_ms2_precursor, const unsigned int& file_counter) { return OpenMS::MetaboTargetedAssay::extractMetaboTargetedAssay(spectra, feature_ms2_index, precursor_rt_tol, precursor_mz_distance, cosine_sim_threshold, transition_threshold, min_fragment_mz, max_fragment_mz, method_consensus_spectrum, exclude_ms2_precursor, file_counter); }, "spectra"_a, "feature_ms2_index"_a, "precursor_rt_tol"_a, "precursor_mz_distance"_a, "cosine_sim_threshold"_a, "transition_threshold"_a, "min_fragment_mz"_a, "max_fragment_mz"_a, "method_consensus_spectrum"_a, "exclude_ms2_precursor"_a, "file_counter"_a, "Extract a vector of MetaboTargetedAssays without using fragment annotation")
        .def_static("extractMetaboTargetedAssayFragmentAnnotation", [](const std::vector<OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair>& v_cmp_spec, const double& transition_threshold, const double& min_fragment_mz, const double& max_fragment_mz, const bool& use_exact_mass, const bool& exclude_ms2_precursor) { return OpenMS::MetaboTargetedAssay::extractMetaboTargetedAssayFragmentAnnotation(v_cmp_spec, transition_threshold, min_fragment_mz, max_fragment_mz, use_exact_mass, exclude_ms2_precursor); }, "v_cmp_spec"_a, "transition_threshold"_a, "min_fragment_mz"_a, "max_fragment_mz"_a, "use_exact_mass"_a, "exclude_ms2_precursor"_a, "Extract a vector of MetaboTargetedAssays using fragment annotation")
        .def_static("pairCompoundWithAnnotatedTDSpectraPairs", [](const std::vector<OpenMS::SiriusMSFile::CompoundInfo>& v_cmpinfo, const std::vector<OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra>& annotated_spectra) { return OpenMS::MetaboTargetedAssay::pairCompoundWithAnnotatedTDSpectraPairs(v_cmpinfo, annotated_spectra); }, "v_cmpinfo"_a, "annotated_spectra"_a, "Pair compound information with annotated target and decoy spectra")
        ;

    // -----------------------------------------------------------------------
    // MetaboTargetedAssay_CompoundTargetDecoyPair
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair>(m, "MetaboTargetedAssay_CompoundTargetDecoyPair", "Compound target-decoy pair for metabolite targeted assays")
        .def(nb::init<>())
        .def(nb::init<OpenMS::SiriusMSFile::CompoundInfo, OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra>())
        .def("__copy__", [](const OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair& self) { return OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair(self); })
        .def("__deepcopy__", [](const OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair& self, nb::dict) { return OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair(self); }, "memo"_a)
        .def_rw("compound_info", &OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair::compound_info)
        .def_rw("target_decoy_spectra", &OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair::target_decoy_spectra)
        ;

    // -----------------------------------------------------------------------
    // MetaboTargetedTargetDecoy
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaboTargetedTargetDecoy>(m, "MetaboTargetedTargetDecoy", 
        R"doc(
Resolve overlapping fragments and missing decoys for experimental
specific decoy generation in targeted/pseudo targeted metabolomics
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaboTargetedTargetDecoy &>())
        .def("__copy__", [](const OpenMS::MetaboTargetedTargetDecoy& self) { return OpenMS::MetaboTargetedTargetDecoy(self); })
        .def("__deepcopy__", [](const OpenMS::MetaboTargetedTargetDecoy& self, nb::dict) { return OpenMS::MetaboTargetedTargetDecoy(self); }, "memo"_a)
        .def_static("constructTargetDecoyMassMapping", [](const OpenMS::TargetedExperiment& t_exp) { return OpenMS::MetaboTargetedTargetDecoy::constructTargetDecoyMassMapping(t_exp); }, "t_exp"_a)
        .def_static("resolveOverlappingTargetDecoyMassesByDecoyMassShift", [](OpenMS::TargetedExperiment& t_exp, std::vector<OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping> mappings, const double& mass_to_add, const double& mz_tol, const OpenMS::String& mz_tol_unit) {
            OpenMS::MetaboTargetedTargetDecoy::resolveOverlappingTargetDecoyMassesByDecoyMassShift(t_exp, mappings, mass_to_add, mz_tol, mz_tol_unit);
            return mappings;
        }, "t_exp"_a, "mappings"_a, "mass_to_add"_a, "mz_tol"_a, "mz_tol_unit"_a,
            R"doc(
Constructs a mass mapping of targets and decoys using the unique m_id identifier
:param t_exp: TransitionExperiment holds compound and transition information used for the mapping
:returns: Updated mappings
)doc")
        .def_static("generateMissingDecoysByMassShift", [](OpenMS::TargetedExperiment& t_exp, std::vector<OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping> mappings, const double& mass_to_add) {
            OpenMS::MetaboTargetedTargetDecoy::generateMissingDecoysByMassShift(t_exp, mappings, mass_to_add);
            return mappings;
        }, "t_exp"_a, "mappings"_a, "mass_to_add"_a,
            R"doc(
Resolves overlapping target and decoy transition masses by adding a specifiable mass (e.g. CH2) to the overlapping decoy fragment
:param t_exp: TransitionExperiment holds compound and transition information
:param mappings: Map of identifier to target and decoy masses
:param mass_to_add: (e.g. CH2)
:param mz_tol: m/z tolerarance for target and decoy transition masses to be considered overlapping
:param mz_tol_unit: m/z tolerance unit
:returns: Updated mappings
)doc")
        ;

    // -----------------------------------------------------------------------
    // MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping>(m, "MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping", ":param mass_to_add: The maximum number of transitions required per assay")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping &>())
        .def("__copy__", [](const OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping& self) { return OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping(self); })
        .def("__deepcopy__", [](const OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping& self, nb::dict) { return OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping(self); }, "memo"_a)
        .def_rw("identifier", &OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping::identifier)
        .def_rw("target_compound_ref", &OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping::target_compound_ref)
        .def_rw("decoy_compound_ref", &OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping::decoy_compound_ref)
        .def_rw("target_product_masses", &OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping::target_product_masses)
        .def_rw("decoy_product_masses", &OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping::decoy_product_masses)
        ;

    // -----------------------------------------------------------------------
    // MorpheusScore
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MorpheusScore>(m, "MorpheusScore", "OpenMS class MorpheusScore")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MorpheusScore &>())
        .def("__copy__", [](const OpenMS::MorpheusScore& self) { return OpenMS::MorpheusScore(self); })
        .def("__deepcopy__", [](const OpenMS::MorpheusScore& self, nb::dict) { return OpenMS::MorpheusScore(self); }, "memo"_a)
        .def_static("compute", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::MSSpectrum& theo_spectrum) { return OpenMS::MorpheusScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum); }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "exp_spectrum"_a, "theo_spectrum"_a, "Returns Morpheus Score")
        .def_static("compute", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::DataArrays::IntegerDataArray& exp_charges, const OpenMS::MSSpectrum& theo_spectrum, const OpenMS::DataArrays::IntegerDataArray& theo_charges) { return OpenMS::MorpheusScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, exp_charges, theo_spectrum, theo_charges); }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "exp_spectrum"_a, "exp_charges"_a, "theo_spectrum"_a, "theo_charges"_a, "Returns Morpheus Score")
        ;

    // -----------------------------------------------------------------------
    // MorpheusScore_Result
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MorpheusScore::Result>(m, "MorpheusScore_Result", "OpenMS class MorpheusScore_Result")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MorpheusScore::Result &>())
        .def("__copy__", [](const OpenMS::MorpheusScore::Result& self) { return OpenMS::MorpheusScore::Result(self); })
        .def("__deepcopy__", [](const OpenMS::MorpheusScore::Result& self, nb::dict) { return OpenMS::MorpheusScore::Result(self); }, "memo"_a)
        .def_rw("matches", &OpenMS::MorpheusScore::Result::matches)
        .def_rw("n_peaks", &OpenMS::MorpheusScore::Result::n_peaks)
        .def_rw("score", &OpenMS::MorpheusScore::Result::score)
        .def_rw("MIC", &OpenMS::MorpheusScore::Result::MIC)
        .def_rw("TIC", &OpenMS::MorpheusScore::Result::TIC)
        .def_rw("err", &OpenMS::MorpheusScore::Result::err)
        .def_rw("err_ppm", &OpenMS::MorpheusScore::Result::err_ppm)
        ;

    // -----------------------------------------------------------------------
    // OPXLDataStructs
    // -----------------------------------------------------------------------
    auto opxldatastructs_class = nb::class_<OpenMS::OPXLDataStructs>(m, "OPXLDataStructs", "OpenMS class OPXLDataStructs")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OPXLDataStructs &>())
        .def("__copy__", [](const OpenMS::OPXLDataStructs& self) { return OpenMS::OPXLDataStructs(self); })
        .def("__deepcopy__", [](const OpenMS::OPXLDataStructs& self, nb::dict) { return OpenMS::OPXLDataStructs(self); }, "memo"_a)
        ;
    // ProteinProteinCrossLinkType enum nested under OPXLDataStructs
    nb::enum_<OpenMS::OPXLDataStructs::ProteinProteinCrossLinkType>(opxldatastructs_class, "ProteinProteinCrossLinkType", nb::is_arithmetic())
        .value("CROSS", OpenMS::OPXLDataStructs::ProteinProteinCrossLinkType::CROSS)
        .value("MONO", OpenMS::OPXLDataStructs::ProteinProteinCrossLinkType::MONO)
        .value("LOOP", OpenMS::OPXLDataStructs::ProteinProteinCrossLinkType::LOOP)
        .value("NUMBER_OF_CROSS_LINK_TYPES", OpenMS::OPXLDataStructs::ProteinProteinCrossLinkType::NUMBER_OF_CROSS_LINK_TYPES)
        .export_values();
    // PeptidePosition enum nested under OPXLDataStructs
    nb::enum_<OpenMS::OPXLDataStructs::PeptidePosition>(opxldatastructs_class, "PeptidePosition", nb::is_arithmetic())
        .value("INTERNAL", OpenMS::OPXLDataStructs::PeptidePosition::INTERNAL)
        .value("C_TERM", OpenMS::OPXLDataStructs::PeptidePosition::C_TERM)
        .value("N_TERM", OpenMS::OPXLDataStructs::PeptidePosition::N_TERM)
        .export_values();

    // --- XLPrecursor ---
    nb::class_<OpenMS::OPXLDataStructs::XLPrecursor>(m, "XLPrecursor", "Cross-link precursor candidate")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::OPXLDataStructs::XLPrecursor& self) { return OpenMS::OPXLDataStructs::XLPrecursor(self); })
        .def("__deepcopy__", [](const OpenMS::OPXLDataStructs::XLPrecursor& self, nb::dict) { return OpenMS::OPXLDataStructs::XLPrecursor(self); }, "memo"_a)
        .def_rw("precursor_mass", &OpenMS::OPXLDataStructs::XLPrecursor::precursor_mass)
        .def_rw("alpha_index", &OpenMS::OPXLDataStructs::XLPrecursor::alpha_index)
        .def_rw("beta_index", &OpenMS::OPXLDataStructs::XLPrecursor::beta_index)
        .def_rw("alpha_seq", &OpenMS::OPXLDataStructs::XLPrecursor::alpha_seq)
        .def_rw("beta_seq", &OpenMS::OPXLDataStructs::XLPrecursor::beta_seq)
        ;

    // --- AASeqWithMass ---
    nb::class_<OpenMS::OPXLDataStructs::AASeqWithMass>(m, "AASeqWithMass", "Amino acid sequence with associated mass")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::OPXLDataStructs::AASeqWithMass& self) { return OpenMS::OPXLDataStructs::AASeqWithMass(self); })
        .def("__deepcopy__", [](const OpenMS::OPXLDataStructs::AASeqWithMass& self, nb::dict) { return OpenMS::OPXLDataStructs::AASeqWithMass(self); }, "memo"_a)
        .def_rw("peptide_mass", &OpenMS::OPXLDataStructs::AASeqWithMass::peptide_mass)
        .def_rw("peptide_seq", &OpenMS::OPXLDataStructs::AASeqWithMass::peptide_seq)
        .def_rw("position", &OpenMS::OPXLDataStructs::AASeqWithMass::position)
        .def_rw("unmodified_seq", &OpenMS::OPXLDataStructs::AASeqWithMass::unmodified_seq)
        ;

    // --- ProteinProteinCrossLink ---
    nb::class_<OpenMS::OPXLDataStructs::ProteinProteinCrossLink>(m, "ProteinProteinCrossLink", "Represents a cross-link between two peptides")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::OPXLDataStructs::ProteinProteinCrossLink& self) { return OpenMS::OPXLDataStructs::ProteinProteinCrossLink(self); })
        .def("__deepcopy__", [](const OpenMS::OPXLDataStructs::ProteinProteinCrossLink& self, nb::dict) { return OpenMS::OPXLDataStructs::ProteinProteinCrossLink(self); }, "memo"_a)
        .def_prop_ro("alpha", [](const OpenMS::OPXLDataStructs::ProteinProteinCrossLink& self) {
            return self.alpha ? *self.alpha : OpenMS::AASequence();
        })
        .def_prop_ro("beta", [](const OpenMS::OPXLDataStructs::ProteinProteinCrossLink& self) {
            return self.beta ? *self.beta : OpenMS::AASequence();
        })
        .def_rw("cross_link_position", &OpenMS::OPXLDataStructs::ProteinProteinCrossLink::cross_link_position)
        .def_rw("cross_linker_mass", &OpenMS::OPXLDataStructs::ProteinProteinCrossLink::cross_linker_mass)
        .def_rw("cross_linker_name", &OpenMS::OPXLDataStructs::ProteinProteinCrossLink::cross_linker_name)
        .def_rw("term_spec_alpha", &OpenMS::OPXLDataStructs::ProteinProteinCrossLink::term_spec_alpha)
        .def_rw("term_spec_beta", &OpenMS::OPXLDataStructs::ProteinProteinCrossLink::term_spec_beta)
        .def_rw("precursor_correction", &OpenMS::OPXLDataStructs::ProteinProteinCrossLink::precursor_correction)
        .def("getType", &OpenMS::OPXLDataStructs::ProteinProteinCrossLink::getType)
        ;

    // --- CrossLinkSpectrumMatch ---
    nb::class_<OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch>(m, "CrossLinkSpectrumMatch", "Result of a cross-link spectrum match")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch& self) { return OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch(self); })
        .def("__deepcopy__", [](const OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch& self, nb::dict) { return OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch(self); }, "memo"_a)
        .def_rw("cross_link", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::cross_link)
        .def_rw("scan_index_light", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::scan_index_light)
        .def_rw("scan_index_heavy", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::scan_index_heavy)
        .def_rw("score", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::score)
        .def_rw("rank", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::rank)
        .def_rw("xquest_score", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::xquest_score)
        .def_rw("pre_score", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::pre_score)
        .def_rw("percTIC", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::percTIC)
        .def_rw("wTIC", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::wTIC)
        .def_rw("wTICold", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::wTICold)
        .def_rw("int_sum", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::int_sum)
        .def_rw("intsum_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::intsum_alpha)
        .def_rw("intsum_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::intsum_beta)
        .def_rw("total_current", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::total_current)
        .def_rw("precursor_error_ppm", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::precursor_error_ppm)
        .def_rw("match_odds", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::match_odds)
        .def_rw("match_odds_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::match_odds_alpha)
        .def_rw("match_odds_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::match_odds_beta)
        .def_rw("log_occupancy", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::log_occupancy)
        .def_rw("log_occupancy_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::log_occupancy_alpha)
        .def_rw("log_occupancy_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::log_occupancy_beta)
        .def_rw("xcorrx_max", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::xcorrx_max)
        .def_rw("xcorrc_max", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::xcorrc_max)
        .def_rw("matched_linear_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::matched_linear_alpha)
        .def_rw("matched_linear_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::matched_linear_beta)
        .def_rw("matched_xlink_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::matched_xlink_alpha)
        .def_rw("matched_xlink_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::matched_xlink_beta)
        .def_rw("num_iso_peaks_mean", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::num_iso_peaks_mean)
        .def_rw("num_iso_peaks_mean_linear_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::num_iso_peaks_mean_linear_alpha)
        .def_rw("num_iso_peaks_mean_linear_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::num_iso_peaks_mean_linear_beta)
        .def_rw("num_iso_peaks_mean_xlinks_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::num_iso_peaks_mean_xlinks_alpha)
        .def_rw("num_iso_peaks_mean_xlinks_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::num_iso_peaks_mean_xlinks_beta)
        .def_rw("ppm_error_abs_sum_linear_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_linear_alpha)
        .def_rw("ppm_error_abs_sum_linear_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_linear_beta)
        .def_rw("ppm_error_abs_sum_xlinks_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_xlinks_alpha)
        .def_rw("ppm_error_abs_sum_xlinks_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_xlinks_beta)
        .def_rw("ppm_error_abs_sum_linear", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_linear)
        .def_rw("ppm_error_abs_sum_xlinks", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_xlinks)
        .def_rw("ppm_error_abs_sum_alpha", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_alpha)
        .def_rw("ppm_error_abs_sum_beta", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum_beta)
        .def_rw("ppm_error_abs_sum", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::ppm_error_abs_sum)
        .def_rw("precursor_correction", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::precursor_correction)
        .def_rw("precursor_total_intensity", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::precursor_total_intensity)
        .def_rw("precursor_target_intensity", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::precursor_target_intensity)
        .def_rw("precursor_signal_proportion", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::precursor_signal_proportion)
        .def_rw("precursor_target_peak_count", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::precursor_target_peak_count)
        .def_rw("precursor_residual_peak_count", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::precursor_residual_peak_count)
        .def_rw("frag_annotations", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::frag_annotations)
        .def_rw("peptide_id_index", &OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch::peptide_id_index)
        ;

    // --- PreprocessedPairSpectra ---
    nb::class_<OpenMS::OPXLDataStructs::PreprocessedPairSpectra>(m, "PreprocessedPairSpectra", "Preprocessed pair spectra for cross-link analysis")
        .def(nb::init<OpenMS::Size>(), "size"_a)
        .def_rw("spectra_linear_peaks", &OpenMS::OPXLDataStructs::PreprocessedPairSpectra::spectra_linear_peaks)
        .def_rw("spectra_xlink_peaks", &OpenMS::OPXLDataStructs::PreprocessedPairSpectra::spectra_xlink_peaks)
        .def_rw("spectra_all_peaks", &OpenMS::OPXLDataStructs::PreprocessedPairSpectra::spectra_all_peaks)
        ;

    // -----------------------------------------------------------------------
    // OPXLHelper
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OPXLHelper>(m, "OPXLHelper", 
        R"doc(
The OPXLHelper class contains functions needed by OpenPepXL to reduce
duplicated code
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OPXLHelper &>())
        .def("__copy__", [](const OpenMS::OPXLHelper& self) { return OpenMS::OPXLHelper(self); })
        .def("__deepcopy__", [](const OpenMS::OPXLHelper& self, nb::dict) { return OpenMS::OPXLHelper(self); }, "memo"_a)
        .def_static("addXLTargetDecoyMV", [](std::vector<OpenMS::PeptideIdentification> peptide_ids) {
            OpenMS::OPXLHelper::addXLTargetDecoyMV(peptide_ids);
            return peptide_ids;
        }, "peptide_ids"_a)
        .def_static("addXLTargetDecoyMV", [](OpenMS::PeptideIdentificationList& peptide_ids) { return OpenMS::OPXLHelper::addXLTargetDecoyMV(peptide_ids); }, "peptide_ids"_a)
        .def_static("addBetaAccessions", [](std::vector<OpenMS::PeptideIdentification> peptide_ids) {
            OpenMS::OPXLHelper::addBetaAccessions(peptide_ids);
            return peptide_ids;
        }, "peptide_ids"_a)
        .def_static("addBetaAccessions", [](OpenMS::PeptideIdentificationList& peptide_ids) { return OpenMS::OPXLHelper::addBetaAccessions(peptide_ids); }, "peptide_ids"_a)
        .def_static("removeBetaPeptideHits", [](std::vector<OpenMS::PeptideIdentification> peptide_ids) {
            OpenMS::OPXLHelper::removeBetaPeptideHits(peptide_ids);
            return peptide_ids;
        }, "peptide_ids"_a)
        .def_static("removeBetaPeptideHits", [](OpenMS::PeptideIdentificationList& peptide_ids) { return OpenMS::OPXLHelper::removeBetaPeptideHits(peptide_ids); }, "peptide_ids"_a)
        .def_static("addPercolatorFeatureList", [](OpenMS::ProteinIdentification& prot_id) { return OpenMS::OPXLHelper::addPercolatorFeatureList(prot_id); }, "prot_id"_a)
        .def_static("computeDeltaScores", []() { std::vector<OpenMS::PeptideIdentification> peptide_ids; OpenMS::OPXLHelper::computeDeltaScores(peptide_ids); return peptide_ids; })
        .def_static("computeDeltaScores", []() { OpenMS::PeptideIdentificationList peptide_ids; OpenMS::OPXLHelper::computeDeltaScores(peptide_ids); return peptide_ids; })
        .def_static("combineTopRanksFromPairs", [](std::vector<OpenMS::PeptideIdentification> peptide_ids, size_t number_top_hits) {
            return OpenMS::OPXLHelper::combineTopRanksFromPairs(peptide_ids, number_top_hits);
        }, "peptide_ids"_a, "number_top_hits"_a)
        .def_static("combineTopRanksFromPairs", [](OpenMS::PeptideIdentificationList& peptide_ids, size_t number_top_hits) { return OpenMS::OPXLHelper::combineTopRanksFromPairs(peptide_ids, number_top_hits); }, "peptide_ids"_a, "number_top_hits"_a)
        .def_static("computePrecursorError", [](const OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch& csm, double precursor_mz, int precursor_charge) { return OpenMS::OPXLHelper::computePrecursorError(csm, precursor_mz, precursor_charge); }, "csm"_a, "precursor_mz"_a, "precursor_charge"_a)
        .def_static("addProteinPositionMetaValues", [](std::vector<OpenMS::PeptideIdentification> peptide_ids) {
            OpenMS::OPXLHelper::addProteinPositionMetaValues(peptide_ids);
            return peptide_ids;
        }, "peptide_ids"_a, "Adds MetaValues for cross-link positions to PeptideHits")
        .def_static("addProteinPositionMetaValues", [](OpenMS::PeptideIdentificationList& peptide_ids) { return OpenMS::OPXLHelper::addProteinPositionMetaValues(peptide_ids); }, "peptide_ids"_a)
        .def_static("isoPeakMeans", [](OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch& csm, const OpenMS::DataArrays::IntegerDataArray& num_iso_peaks_array, const std::vector<std::pair<size_t, size_t>>& matched_spec_linear_alpha, const std::vector<std::pair<size_t, size_t>>& matched_spec_linear_beta, const std::vector<std::pair<size_t, size_t>>& matched_spec_xlinks_alpha, const std::vector<std::pair<size_t, size_t>>& matched_spec_xlinks_beta) { OpenMS::OPXLHelper::isoPeakMeans(csm, num_iso_peaks_array, matched_spec_linear_alpha, matched_spec_linear_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta); }, "csm"_a, "num_iso_peaks_array"_a, "matched_spec_linear_alpha"_a, "matched_spec_linear_beta"_a, "matched_spec_xlinks_alpha"_a, "matched_spec_xlinks_beta"_a, "Computes the mean of alpha, beta, xlinks-alpha and xlinks-beta iso peak counts")
        .def_static("enumerateCrossLinksAndMasses", [](
                const std::vector<OpenMS::OPXLDataStructs::AASeqWithMass>& peptides,
                double cross_link_mass_light,
                const OpenMS::DoubleList& cross_link_mass_mono_link,
                const OpenMS::StringList& cross_link_residue1,
                const OpenMS::StringList& cross_link_residue2,
                const std::vector<double>& spectrum_precursors,
                double precursor_mass_tolerance,
                bool precursor_mass_tolerance_unit_ppm) {
            std::vector<int> precursor_correction_positions;
            auto result = OpenMS::OPXLHelper::enumerateCrossLinksAndMasses(
                peptides, cross_link_mass_light, cross_link_mass_mono_link,
                cross_link_residue1, cross_link_residue2, spectrum_precursors,
                precursor_correction_positions, precursor_mass_tolerance,
                precursor_mass_tolerance_unit_ppm);
            return nb::make_tuple(result, precursor_correction_positions);
        }, "peptides"_a, "cross_link_mass_light"_a, "cross_link_mass_mono_link"_a,
           "cross_link_residue1"_a, "cross_link_residue2"_a, "spectrum_precursors"_a,
           "precursor_mass_tolerance"_a, "precursor_mass_tolerance_unit_ppm"_a,
           "Enumerates cross-link candidates and masses, returns (candidates, correction_positions)")
        .def_static("digestDatabase", [](
                std::vector<OpenMS::FASTAFile::FASTAEntry> fasta_db,
                const OpenMS::EnzymaticDigestion& digestor,
                size_t min_peptide_length,
                const OpenMS::StringList& cross_link_residue1,
                const OpenMS::StringList& cross_link_residue2,
                const OpenMS::ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
                const OpenMS::ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
                size_t max_variable_mods_per_peptide) {
            return OpenMS::OPXLHelper::digestDatabase(fasta_db, digestor, min_peptide_length,
                cross_link_residue1, cross_link_residue2, fixed_modifications,
                variable_modifications, max_variable_mods_per_peptide);
        }, "fasta_db"_a, "digestor"_a, "min_peptide_length"_a,
           "cross_link_residue1"_a, "cross_link_residue2"_a,
           "fixed_modifications"_a, "variable_modifications"_a,
           "max_variable_mods_per_peptide"_a,
           "Digests a FASTA database and returns peptides with masses")
        .def_static("buildCandidates", [](
                const std::vector<OpenMS::OPXLDataStructs::XLPrecursor>& candidates,
                const std::vector<int>& precursor_corrections,
                const std::vector<int>& precursor_correction_positions,
                const std::vector<OpenMS::OPXLDataStructs::AASeqWithMass>& peptide_masses,
                const OpenMS::StringList& cross_link_residue1,
                const OpenMS::StringList& cross_link_residue2,
                double cross_link_mass,
                const OpenMS::DoubleList& cross_link_mass_mono_link,
                const std::vector<double>& spectrum_precursor_vector,
                const std::vector<double>& allowed_error_vector,
                const OpenMS::String& cross_link_name) {
            return OpenMS::OPXLHelper::buildCandidates(candidates, precursor_corrections,
                precursor_correction_positions, peptide_masses, cross_link_residue1,
                cross_link_residue2, cross_link_mass, cross_link_mass_mono_link,
                spectrum_precursor_vector, allowed_error_vector, cross_link_name);
        }, "candidates"_a, "precursor_corrections"_a, "precursor_correction_positions"_a,
           "peptide_masses"_a, "cross_link_residue1"_a, "cross_link_residue2"_a,
           "cross_link_mass"_a, "cross_link_mass_mono_link"_a,
           "spectrum_precursor_vector"_a, "allowed_error_vector"_a, "cross_link_name"_a,
           "Builds cross-link candidates from precursor candidates")
        .def_static("buildFragmentAnnotations", [](
                const std::vector<std::pair<size_t, size_t>>& matching,
                const OpenMS::MSSpectrum& theoretical_spectrum,
                const OpenMS::MSSpectrum& experiment_spectrum) {
            std::vector<OpenMS::PeptideHit::PeakAnnotation> frag_annotations;
            OpenMS::OPXLHelper::buildFragmentAnnotations(frag_annotations, matching,
                theoretical_spectrum, experiment_spectrum);
            return frag_annotations;
        }, "matching"_a, "theoretical_spectrum"_a, "experiment_spectrum"_a,
           "Builds fragment annotations from spectrum matching results")
        .def_static("buildPeptideIDs", [](
                std::vector<OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch>& top_csms_spectrum,
                std::vector<std::vector<OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch>>& all_top_csms,
                size_t all_top_csms_current_index,
                const OpenMS::PeakMap& spectra,
                size_t scan_index,
                size_t scan_index_heavy) {
            OpenMS::PeptideIdentificationList peptide_ids;
            OpenMS::OPXLHelper::buildPeptideIDs(peptide_ids, top_csms_spectrum, all_top_csms,
                all_top_csms_current_index, spectra, scan_index, scan_index_heavy);
            return peptide_ids;
        }, "top_csms_spectrum"_a, "all_top_csms"_a, "all_top_csms_current_index"_a,
           "spectra"_a, "scan_index"_a, "scan_index_heavy"_a,
           "Builds PeptideIdentifications from cross-link spectrum matches")
        .def_static("collectPrecursorCandidates", [](
                const OpenMS::IntList& precursor_correction_steps,
                double precursor_mass,
                double precursor_mass_tolerance,
                bool precursor_mass_tolerance_unit_ppm,
                const std::vector<OpenMS::OPXLDataStructs::AASeqWithMass>& filtered_peptide_masses,
                double cross_link_mass,
                const OpenMS::DoubleList& cross_link_mass_mono_link,
                const OpenMS::StringList& cross_link_residue1,
                const OpenMS::StringList& cross_link_residue2,
                const OpenMS::String& cross_link_name,
                bool use_sequence_tags,
                const std::vector<std::string>& tags) {
            return OpenMS::OPXLHelper::collectPrecursorCandidates(
                precursor_correction_steps, precursor_mass, precursor_mass_tolerance,
                precursor_mass_tolerance_unit_ppm, filtered_peptide_masses, cross_link_mass,
                cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2,
                cross_link_name, use_sequence_tags, tags);
        }, "precursor_correction_steps"_a, "precursor_mass"_a,
           "precursor_mass_tolerance"_a, "precursor_mass_tolerance_unit_ppm"_a,
           "filtered_peptide_masses"_a, "cross_link_mass"_a,
           "cross_link_mass_mono_link"_a, "cross_link_residue1"_a,
           "cross_link_residue2"_a, "cross_link_name"_a,
           "use_sequence_tags"_a = false, "tags"_a = std::vector<std::string>(),
           "Collects precursor candidates for a given precursor mass")
        ;

    // -----------------------------------------------------------------------
    // OPXLSpectrumProcessingAlgorithms
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OPXLSpectrumProcessingAlgorithms>(m, "OPXLSpectrumProcessingAlgorithms", "OpenMS class OPXLSpectrumProcessingAlgorithms")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OPXLSpectrumProcessingAlgorithms &>())
        .def("__copy__", [](const OpenMS::OPXLSpectrumProcessingAlgorithms& self) { return OpenMS::OPXLSpectrumProcessingAlgorithms(self); })
        .def("__deepcopy__", [](const OpenMS::OPXLSpectrumProcessingAlgorithms& self, nb::dict) { return OpenMS::OPXLSpectrumProcessingAlgorithms(self); }, "memo"_a)
        .def_static("mergeAnnotatedSpectra", [](OpenMS::MSSpectrum& first_spectrum, OpenMS::MSSpectrum& second_spectrum) { return OpenMS::OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(first_spectrum, second_spectrum); }, "first_spectrum"_a, "second_spectrum"_a)
        .def_static("preprocessSpectra", [](OpenMS::MSExperiment& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, size_t peptide_min_size, int min_precursor_charge, int max_precursor_charge, bool deisotope, bool labeled) { return OpenMS::OPXLSpectrumProcessingAlgorithms::preprocessSpectra(exp, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peptide_min_size, min_precursor_charge, max_precursor_charge, deisotope, labeled); }, "exp"_a, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "peptide_min_size"_a, "min_precursor_charge"_a, "max_precursor_charge"_a, "deisotope"_a, "labeled"_a)
        .def_static("getSpectrumAlignmentFastCharge", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const OpenMS::MSSpectrum& theo_spectrum, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::DataArrays::IntegerDataArray& theo_charges, const OpenMS::DataArrays::IntegerDataArray& exp_charges, double intensity_cutoff) {
            std::vector<std::pair<size_t, size_t>> alignment;
            OpenMS::DataArrays::FloatDataArray ppm_error_array;
            OpenMS::OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, theo_spectrum, exp_spectrum, theo_charges, exp_charges, ppm_error_array, intensity_cutoff);
            return std::make_pair(alignment, ppm_error_array);
        }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "theo_spectrum"_a, "exp_spectrum"_a, "theo_charges"_a, "exp_charges"_a, "intensity_cutoff"_a = 0.0, "Computes a spectrum alignment considering fragment charges, returns (alignment, ppm_error_array)")
        .def_static("getSpectrumAlignmentSimple", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::vector<OpenMS::SimpleTSGXLMS::SimplePeak>& theo_spectrum, const OpenMS::MSSpectrum& exp_spectrum, const OpenMS::DataArrays::IntegerDataArray& exp_charges) {
            std::vector<std::pair<size_t, size_t>> alignment;
            OpenMS::OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(alignment, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, theo_spectrum, exp_spectrum, exp_charges);
            return alignment;
        }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "theo_spectrum"_a, "exp_spectrum"_a, "exp_charges"_a, "Computes a spectrum alignment using SimplePeak for the theoretical spectrum")
        ;

    // -----------------------------------------------------------------------
    // OpenSwathDataAccessHelper
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSwathDataAccessHelper>(m, "OpenSwathDataAccessHelper", 
        R"doc(
Several helpers to convert OpenMS datastructures to structures that
implement the OpenSWATH interfaces
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSwathDataAccessHelper &>())
        .def("__copy__", [](const OpenMS::OpenSwathDataAccessHelper& self) { return OpenMS::OpenSwathDataAccessHelper(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSwathDataAccessHelper& self, nb::dict) { return OpenMS::OpenSwathDataAccessHelper(self); }, "memo"_a)
        .def_static("convertToOpenMSSpectrum", [](const std::shared_ptr<OpenSwath::OSSpectrum>& sptr, OpenMS::MSSpectrum& spectrum) { return OpenMS::OpenSwathDataAccessHelper::convertToOpenMSSpectrum(sptr, spectrum); }, "sptr"_a, "spectrum"_a, "Converts a SpectrumPtr to an OpenMS Spectrum")
        .def_static("convertToOpenMSChromatogram", [](const std::shared_ptr<OpenSwath::OSChromatogram>& cptr, OpenMS::MSChromatogram& chromatogram) { return OpenMS::OpenSwathDataAccessHelper::convertToOpenMSChromatogram(cptr, chromatogram); }, "cptr"_a, "chromatogram"_a, "Converts a ChromatogramPtr to an OpenMS Chromatogram")
        .def_static("convertToOpenMSChromatogramFilter", [](OpenMS::MSChromatogram& chromatogram, const std::shared_ptr<OpenSwath::OSChromatogram>& cptr, double rt_min, double rt_max) { return OpenMS::OpenSwathDataAccessHelper::convertToOpenMSChromatogramFilter(chromatogram, cptr, rt_min, rt_max); }, "chromatogram"_a, "cptr"_a, "rt_min"_a, "rt_max"_a)
        .def_static("convertTargetedExp", [](const OpenMS::TargetedExperiment& transition_exp_, OpenSwath::LightTargetedExperiment& transition_exp) { return OpenMS::OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp); }, "transition_exp_"_a, "transition_exp"_a, "Converts from the OpenMS TargetedExperiment to the OpenMs LightTargetedExperiment")
        .def_static("convertPeptideToAASequence", [](const OpenSwath::LightCompound& peptide, OpenMS::AASequence& aa_sequence) { return OpenMS::OpenSwathDataAccessHelper::convertPeptideToAASequence(peptide, aa_sequence); }, "peptide"_a, "aa_sequence"_a, "Converts from the LightCompound to an OpenMS AASequence (with correct modifications)")
        .def_static("convertTargetedCompound", [](const OpenMS::TargetedExperiment::Peptide& pep) {
            OpenSwath::LightCompound comp;
            OpenMS::OpenSwathDataAccessHelper::convertTargetedCompound(pep, comp);
            return comp;
        }, "pep"_a, "Converts a TargetedExperiment Peptide to a LightCompound")
        .def_static("convertToSpectrumPtr", [](const OpenMS::MSSpectrum& spectrum) { return OpenMS::OpenSwathDataAccessHelper::convertToSpectrumPtr(spectrum); }, "spectrum"_a, "Converts an OpenMS Spectrum to a SpectrumPtr")
        .def_static("convertToChromatogramPtr", [](const OpenMS::MSChromatogram& chromatogram) { return OpenMS::OpenSwathDataAccessHelper::convertToChromatogramPtr(chromatogram); }, "chromatogram"_a, "Converts an OpenMS Chromatogram to a ChromatogramPtr")
        ;

    // -----------------------------------------------------------------------
    // OpenSwathHelper
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSwathHelper>(m, "OpenSwathHelper", "A helper class that is used by several OpenSWATH tools")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OpenSwathHelper &>())
        .def("__copy__", [](const OpenMS::OpenSwathHelper& self) { return OpenMS::OpenSwathHelper(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSwathHelper& self, nb::dict) { return OpenMS::OpenSwathHelper(self); }, "memo"_a)
        .def_static("computePrecursorId", [](const OpenMS::String& transition_group_id, int isotope) { return OpenMS::OpenSwathHelper::computePrecursorId(transition_group_id, isotope); }, "transition_group_id"_a, "isotope"_a, 
            R"doc(
Computes the min and max retention time value
Estimate the retention time span of a targeted experiment by returning the min/max values in retention time as a pair
:return: A std `pair` that contains (min,max)
)doc")
        .def_static("estimateRTRange", [](const OpenSwath::LightTargetedExperiment& exp) { return OpenMS::OpenSwathHelper::estimateRTRange(exp); }, "exp"_a)
        .def_static("checkSwathMapAndSelectTransitions", [](const OpenMS::PeakMap& exp, const OpenMS::TargetedExperiment& targeted_exp, double min_upper_edge_dist) {
            OpenMS::TargetedExperiment selected_transitions;
            bool result = OpenMS::OpenSwathHelper::checkSwathMapAndSelectTransitions(exp, targeted_exp, selected_transitions, min_upper_edge_dist);
            return nb::make_tuple(result, selected_transitions);
        }, "exp"_a, "targeted_exp"_a, "min_upper_edge_dist"_a, "Checks a swath map and selects appropriate transitions. Returns (success, selected_transitions)")
        ;

    // -----------------------------------------------------------------------
    // OpenSwathOSWWriter
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSwathOSWWriter>(m, "OpenSwathOSWWriter", "Class to write out an OpenSwath OSW SQLite output (PyProphet input)")
        .def(nb::init<OpenMS::String, bool>())
        .def("__copy__", [](const OpenMS::OpenSwathOSWWriter& self) { return OpenMS::OpenSwathOSWWriter(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSwathOSWWriter& self, nb::dict) { return OpenMS::OpenSwathOSWWriter(self); }, "memo"_a)
        .def("isActive", [](const OpenMS::OpenSwathOSWWriter& self) { return self.isActive(); })
        .def("writeHeader", [](OpenMS::OpenSwathOSWWriter& self) { return self.writeHeader(); }, "Initializes file by generating SQLite tables")
        .def("addRun", [](OpenMS::OpenSwathOSWWriter& self, size_t run_id, const OpenMS::String& input_filename) { return self.addRun(run_id, input_filename); }, "run_id"_a, "input_filename"_a, 
            R"doc(
Write data to disk
Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
:param to_osw_output: Statements generated by prepareLine
)doc")
        .def("setRunId", [](OpenMS::OpenSwathOSWWriter& self, size_t run_id) { return self.setRunId(run_id); }, "run_id"_a)
        .def("writeLines", [](OpenMS::OpenSwathOSWWriter& self, const std::vector<OpenMS::String>& to_osw_output) { return self.writeLines(to_osw_output); }, "to_osw_output"_a, 
            R"doc(
Prepare a single line (feature) for output
The result can be flushed to disk using writeLines (either line by line or after collecting several lines)
:param pep: The compound (peptide/metabolite) used for extraction
:param transition: The transition used for extraction
:param output: The feature map containing all features (each feature will generate one entry in the output)
:param id: The transition group identifier (peptide/metabolite id)
:return: A String to be written using writeLines
)doc")
        .def("prepareLine", [](const OpenMS::OpenSwathOSWWriter& self, const OpenSwath::LightCompound& pep, const OpenSwath::LightTransition& transition, const OpenMS::FeatureMap& output, const OpenMS::String& id) {
            return self.prepareLine(pep, &transition, output, id);
        }, "pep"_a, "transition"_a, "output"_a, "id"_a, "Prepare a single line (feature) for output")
        ;

    // -----------------------------------------------------------------------
    // PI_PeakArea
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakIntegrator::PeakArea>(m, "PI_PeakArea", "OpenMS class PI_PeakArea")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeakIntegrator::PeakArea &>())
        .def("__copy__", [](const OpenMS::PeakIntegrator::PeakArea& self) { return OpenMS::PeakIntegrator::PeakArea(self); })
        .def("__deepcopy__", [](const OpenMS::PeakIntegrator::PeakArea& self, nb::dict) { return OpenMS::PeakIntegrator::PeakArea(self); }, "memo"_a)
        .def_rw("area", &OpenMS::PeakIntegrator::PeakArea::area)
        .def_rw("height", &OpenMS::PeakIntegrator::PeakArea::height)
        .def_rw("apex_pos", &OpenMS::PeakIntegrator::PeakArea::apex_pos)
        .def_rw("hull_points", &OpenMS::PeakIntegrator::PeakArea::hull_points)
        ;

    // -----------------------------------------------------------------------
    // PI_PeakBackground
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakIntegrator::PeakBackground>(m, "PI_PeakBackground", "OpenMS class PI_PeakBackground")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeakIntegrator::PeakBackground &>())
        .def("__copy__", [](const OpenMS::PeakIntegrator::PeakBackground& self) { return OpenMS::PeakIntegrator::PeakBackground(self); })
        .def("__deepcopy__", [](const OpenMS::PeakIntegrator::PeakBackground& self, nb::dict) { return OpenMS::PeakIntegrator::PeakBackground(self); }, "memo"_a)
        .def_rw("area", &OpenMS::PeakIntegrator::PeakBackground::area)
        .def_rw("height", &OpenMS::PeakIntegrator::PeakBackground::height)
        ;

    // -----------------------------------------------------------------------
    // PI_PeakShapeMetrics
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakIntegrator::PeakShapeMetrics>(m, "PI_PeakShapeMetrics", "OpenMS class PI_PeakShapeMetrics")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeakIntegrator::PeakShapeMetrics &>())
        .def("__copy__", [](const OpenMS::PeakIntegrator::PeakShapeMetrics& self) { return OpenMS::PeakIntegrator::PeakShapeMetrics(self); })
        .def("__deepcopy__", [](const OpenMS::PeakIntegrator::PeakShapeMetrics& self, nb::dict) { return OpenMS::PeakIntegrator::PeakShapeMetrics(self); }, "memo"_a)
        .def_rw("width_at_5", &OpenMS::PeakIntegrator::PeakShapeMetrics::width_at_5)
        .def_rw("width_at_10", &OpenMS::PeakIntegrator::PeakShapeMetrics::width_at_10)
        .def_rw("width_at_50", &OpenMS::PeakIntegrator::PeakShapeMetrics::width_at_50)
        .def_rw("start_position_at_5", &OpenMS::PeakIntegrator::PeakShapeMetrics::start_position_at_5)
        .def_rw("start_position_at_10", &OpenMS::PeakIntegrator::PeakShapeMetrics::start_position_at_10)
        .def_rw("start_position_at_50", &OpenMS::PeakIntegrator::PeakShapeMetrics::start_position_at_50)
        .def_rw("end_position_at_5", &OpenMS::PeakIntegrator::PeakShapeMetrics::end_position_at_5)
        .def_rw("end_position_at_10", &OpenMS::PeakIntegrator::PeakShapeMetrics::end_position_at_10)
        .def_rw("end_position_at_50", &OpenMS::PeakIntegrator::PeakShapeMetrics::end_position_at_50)
        .def_rw("total_width", &OpenMS::PeakIntegrator::PeakShapeMetrics::total_width)
        .def_rw("tailing_factor", &OpenMS::PeakIntegrator::PeakShapeMetrics::tailing_factor)
        .def_rw("asymmetry_factor", &OpenMS::PeakIntegrator::PeakShapeMetrics::asymmetry_factor)
        .def_rw("slope_of_baseline", &OpenMS::PeakIntegrator::PeakShapeMetrics::slope_of_baseline)
        .def_rw("baseline_delta_2_height", &OpenMS::PeakIntegrator::PeakShapeMetrics::baseline_delta_2_height)
        .def_rw("points_across_baseline", &OpenMS::PeakIntegrator::PeakShapeMetrics::points_across_baseline)
        .def_rw("points_across_half_height", &OpenMS::PeakIntegrator::PeakShapeMetrics::points_across_half_height)
        ;

    // -----------------------------------------------------------------------
    // PScore
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PScore>(m, "PScore", "OpenMS class PScore")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PScore &>())
        .def("__copy__", [](const OpenMS::PScore& self) { return OpenMS::PScore(self); })
        .def("__deepcopy__", [](const OpenMS::PScore& self, nb::dict) { return OpenMS::PScore(self); }, "memo"_a)
        .def_static("calculateIntensityRankInMZWindow", [](const std::vector<double>& mz, const std::vector<double>& intensities, double mz_window) { return OpenMS::PScore::calculateIntensityRankInMZWindow(mz, intensities, mz_window); }, "mz"_a, "intensities"_a, "mz_window"_a)
        .def_static("calculateRankMap", [](const OpenMS::MSExperiment& peak_map, double mz_window) { return OpenMS::PScore::calculateRankMap(peak_map, mz_window); }, "peak_map"_a, "mz_window"_a, 
            R"doc(
Calculate local (windowed) peak ranks
The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity
The result can be used to efficiently filter spectra for top 1..n peaks in mass windows
:param mz: The m/z positions of the peaks
:param intensities: The intensities of the peaks
:param mz_window: The window in Thomson centered at each peak
)doc")
        .def_static("calculatePeakLevelSpectra", [](const OpenMS::MSSpectrum& spec, const std::vector<size_t>& ranks, size_t min_level, size_t max_level) { return OpenMS::PScore::calculatePeakLevelSpectra(spec, ranks, min_level, max_level); }, "spec"_a, "ranks"_a, "min_level"_a, "max_level"_a,
            R"doc(
Precalculated, windowed peak ranks for a whole experiment
The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity
:param peak_map: Fragment spectra used for rank calculation. Typically a peak map after removal of all MS1 spectra
:param mz_window: Window in Thomson centered at each peak
)doc")
        .def_static("computePScore", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<size_t, OpenMS::MSSpectrum>& peak_level_spectra, const std::vector<OpenMS::MSSpectrum>& theo_spectra, double mz_window) { return OpenMS::PScore::computePScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peak_level_spectra, theo_spectra, mz_window); }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "peak_level_spectra"_a, "theo_spectra"_a, "mz_window"_a = 100.0, "Compute the PScore for a vector of theoretical spectra")
        .def_static("computePScore", [](double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<size_t, OpenMS::MSSpectrum>& peak_level_spectra, const OpenMS::MSSpectrum& theo_spectrum, double mz_window) { return OpenMS::PScore::computePScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peak_level_spectra, theo_spectrum, mz_window); }, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "peak_level_spectra"_a, "theo_spectrum"_a, "mz_window"_a = 100.0, "Compute the PScore for a single theoretical spectrum")
        ;

    // -----------------------------------------------------------------------
    // PeakGroup
    // -----------------------------------------------------------------------
    auto peakgroup_class = nb::class_<OpenMS::PeakGroup>(m, "PeakGroup", 
        R"doc(
Class describing a deconvolved mass.
A mass contains multiple (LogMz) peaks of different charges and isotope indices.
PeakGroup is the set of such peaks representing a single monoisotopic mass.
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<int, int, bool>())
        .def(nb::init<const OpenMS::PeakGroup &>())
        .def("__copy__", [](const OpenMS::PeakGroup& self) { return OpenMS::PeakGroup(self); })
        .def("__deepcopy__", [](const OpenMS::PeakGroup& self, nb::dict) { return OpenMS::PeakGroup(self); }, "memo"_a)
        .def(nb::self < nb::self)
        .def(nb::self > nb::self)
        .def(nb::self == nb::self)
        .def("setScanNumber", [](OpenMS::PeakGroup& self, int scan_number) { return self.setScanNumber(scan_number); }, "scan_number"_a, "Sets the scan number")
        .def("setChargeIsotopeCosine", [](OpenMS::PeakGroup& self, int abs_charge, float cos) { return self.setChargeIsotopeCosine(abs_charge, cos); }, "abs_charge"_a, "cos"_a, "Sets isotope cosine for given charge")
        .def("setIsotopeCosine", [](OpenMS::PeakGroup& self, float cos) { return self.setIsotopeCosine(cos); }, "cos"_a, "Sets the isotope cosine score")
        .def("setRepAbsCharge", [](OpenMS::PeakGroup& self, int max_snr_abs_charge) { return self.setRepAbsCharge(max_snr_abs_charge); }, "max_snr_abs_charge"_a, "Sets the representative charge")
        .def("setMonoisotopicMass", [](OpenMS::PeakGroup& self, double mono_mass) { return self.setMonoisotopicMass(mono_mass); }, "mono_mass"_a, "Sets the monoisotopic mass")
        .def("setQscore", [](OpenMS::PeakGroup& self, double qscore) { return self.setQscore(qscore); }, "qscore"_a, "Sets the quality score")
        .def("setChargeScore", [](OpenMS::PeakGroup& self, float charge_score) { return self.setChargeScore(charge_score); }, "charge_score"_a, "Sets the charge fit score")
        .def("setAvgPPMError", [](OpenMS::PeakGroup& self, float error) { return self.setAvgPPMError(error); }, "error"_a, "Sets average ppm error")
        .def("setSNR", [](OpenMS::PeakGroup& self, float snr) { return self.setSNR(snr); }, "snr"_a, "Sets the SNR")
        .def("setChargeSNR", [](OpenMS::PeakGroup& self, int abs_charge, float c_snr) { return self.setChargeSNR(abs_charge, c_snr); }, "abs_charge"_a, "c_snr"_a, "Sets SNR for given charge")
        .def("setTargeted", [](OpenMS::PeakGroup& self) { return self.setTargeted(); }, "Marks this peak group as targeted")
        .def("getScanNumber", [](const OpenMS::PeakGroup& self) { return self.getScanNumber(); }, "Returns the scan number")
        .def("getIntensity", [](const OpenMS::PeakGroup& self) { return self.getIntensity(); }, "Returns the summed intensity")
        .def("getChargeSNR", [](const OpenMS::PeakGroup& self, int abs_charge) { return self.getChargeSNR(abs_charge); }, "abs_charge"_a, "Returns SNR for given charge")
        .def("getChargeIsotopeCosine", [](const OpenMS::PeakGroup& self, int abs_charge) { return self.getChargeIsotopeCosine(abs_charge); }, "abs_charge"_a, "Returns isotope cosine for given charge")
        .def("getChargeIntensity", [](const OpenMS::PeakGroup& self, int abs_charge) { return self.getChargeIntensity(abs_charge); }, "abs_charge"_a, "Returns intensity for given charge")
        .def("getIsotopeIntensities", [](const OpenMS::PeakGroup& self) -> const std::vector<float> & { return self.getIsotopeIntensities(); }, nb::rv_policy::reference_internal, "Returns per-isotope intensities")
        .def("getIsotopeCosine", [](const OpenMS::PeakGroup& self) { return self.getIsotopeCosine(); }, "Returns the isotope cosine score")
        .def("getPeakOccupancy", [](const OpenMS::PeakGroup& self) { return self.getPeakOccupancy(); }, "Returns peak occupancy (0-1)")
        .def("getRepAbsCharge", [](const OpenMS::PeakGroup& self) { return self.getRepAbsCharge(); }, "Returns the representative charge")
        .def("getQscore", [](const OpenMS::PeakGroup& self) { return self.getQscore(); }, "Returns the quality score (0-1)")
        .def("getQscore2D", [](const OpenMS::PeakGroup& self) { return self.getQscore2D(); }, "Returns the 2D quality score incorporating feature-level information")
        .def("getSNR", [](const OpenMS::PeakGroup& self) { return self.getSNR(); }, "Returns the signal-to-noise ratio")
        .def("getChargeScore", [](const OpenMS::PeakGroup& self) { return self.getChargeScore(); }, "Returns the charge fit score")
        .def("getAvgPPMError", [](const OpenMS::PeakGroup& self) { return self.getAvgPPMError(); }, "Returns average ppm error")
        .def("getAvgDaError", [](const OpenMS::PeakGroup& self) { return self.getAvgDaError(); }, "Returns average Da error")
        .def("isPositive", [](const OpenMS::PeakGroup& self) { return self.isPositive(); }, "Returns true if positive ionization mode")
        .def("isTargeted", [](const OpenMS::PeakGroup& self) { return self.isTargeted(); }, "Returns true if this peak group was targeted")
        .def("getTargetDecoyType", [](const OpenMS::PeakGroup& self) { return self.getTargetDecoyType(); }, "Returns target/decoy type")
        .def("setTargetDecoyType", [](OpenMS::PeakGroup& self, OpenMS::PeakGroup::TargetDecoyType index) { return self.setTargetDecoyType(index); }, "index"_a, "Sets the target/decoy type")
        .def("getQvalue", [](const OpenMS::PeakGroup& self) { return self.getQvalue(); }, "Returns the q-value for FDR")
        .def("setQvalue", [](OpenMS::PeakGroup& self, double q) { return self.setQvalue(q); }, "q"_a, "Sets the q-value")
        .def("setIsotopeDaDistance", [](OpenMS::PeakGroup& self, double d) { return self.setIsotopeDaDistance(d); }, "d"_a, "Sets isotope distance")
        .def("getIsotopeDaDistance", [](const OpenMS::PeakGroup& self) { return self.getIsotopeDaDistance(); }, "Returns distance between consecutive isotopes")
        .def("getMinNegativeIsotopeIndex", [](const OpenMS::PeakGroup& self) { return self.getMinNegativeIsotopeIndex(); }, "Returns minimum negative isotope index")
        .def("setIndex", [](OpenMS::PeakGroup& self, unsigned int i) { return self.setIndex(i); }, "i"_a, "Sets the peak group index")
        .def("setQscore2D", [](OpenMS::PeakGroup& self, double fqscore) { return self.setQscore2D(fqscore); }, "fqscore"_a, "Sets the 2D quality score")
        .def("setFeatureIndex", [](OpenMS::PeakGroup& self, unsigned int findex) { return self.setFeatureIndex(findex); }, "findex"_a, "Sets the feature index")
        .def("getIndex", [](const OpenMS::PeakGroup& self) { return self.getIndex(); }, "Returns the peak group index")
        .def("getFeatureIndex", [](const OpenMS::PeakGroup& self) { return self.getFeatureIndex(); }, "Returns the feature index")
        .def("begin", [](const OpenMS::PeakGroup& self) { return self.begin(); })
        .def("end", [](const OpenMS::PeakGroup& self) { return self.end(); })
        .def("__getitem__", [](OpenMS::PeakGroup& self, size_t i) { if (i >= self.size()) throw nb::index_error(); return self[i]; })
        .def("getMassErrors", [](const OpenMS::PeakGroup& self, bool ppm) { return self.getMassErrors(ppm); }, "ppm"_a = true, "Returns mass errors per isotope")
        .def("push_back", [](OpenMS::PeakGroup& self, const OpenMS::FLASHHelperClasses::LogMzPeak& pg) { return self.push_back(pg); }, "pg"_a, "Adds a LogMzPeak")
        .def("size", [](const OpenMS::PeakGroup& self) { return self.size(); }, "Returns number of LogMzPeaks")
        .def("reserve", [](OpenMS::PeakGroup& self, size_t n) { return self.reserve(n); }, "n"_a, "Reserves space for n peaks")
        .def("empty", [](const OpenMS::PeakGroup& self) { return self.empty(); }, "Returns true if no peaks")
        .def("sort", [](OpenMS::PeakGroup& self) { return self.sort(); }, "Sorts peaks by log m/z")
        .def("__iter__", [](OpenMS::PeakGroup& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::PeakGroup>(), "PeakGroup_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::PeakGroup& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::PeakGroup& self, size_t i) -> const OpenMS::FLASHHelperClasses::LogMzPeak & { 
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)

        .def("getMonoMass", &OpenMS::PeakGroup::getMonoMass, "Returns the monoisotopic mass")
        ;
    // TargetDecoyType enum nested under PeakGroup
    nb::enum_<OpenMS::PeakGroup::TargetDecoyType>(peakgroup_class, "TargetDecoyType", nb::is_arithmetic())
        .value("target", OpenMS::PeakGroup::TargetDecoyType::target)
        .value("noise_decoy", OpenMS::PeakGroup::TargetDecoyType::noise_decoy)
        .value("signal_decoy", OpenMS::PeakGroup::TargetDecoyType::signal_decoy)
        .export_values();

    // -----------------------------------------------------------------------
    // TargetedExperiment_Modification (TargetedExperimentHelper::Peptide::Modification)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::Peptide::Modification, OpenMS::CVTermListInterface>(m, "TargetedExperiment_Modification",
        "Modification on a peptide in a targeted experiment")
        .def(nb::init<>())
        .def_rw("avg_mass_delta", &OpenMS::TargetedExperimentHelper::Peptide::Modification::avg_mass_delta)
        .def_rw("mono_mass_delta", &OpenMS::TargetedExperimentHelper::Peptide::Modification::mono_mass_delta)
        .def_rw("location", &OpenMS::TargetedExperimentHelper::Peptide::Modification::location)
        .def_rw("unimod_id", &OpenMS::TargetedExperimentHelper::Peptide::Modification::unimod_id)
        ;

    // -----------------------------------------------------------------------
    // Peptide
    // -----------------------------------------------------------------------
    auto peptide_class = nb::class_<OpenMS::TargetedExperimentHelper::Peptide>(m, "Peptide", 
        R"doc(
Represents a peptide in a targeted experiment (e.g. used by
ReactionMonitoringTransition and IncludeExcludeTarget)
CVTermList
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TargetedExperimentHelper::Peptide &>())
        .def("__copy__", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return OpenMS::TargetedExperimentHelper::Peptide(self); })
        .def("__deepcopy__", [](const OpenMS::TargetedExperimentHelper::Peptide& self, nb::dict) { return OpenMS::TargetedExperimentHelper::Peptide(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("setPeptideGroupLabel", [](OpenMS::TargetedExperimentHelper::Peptide& self, const OpenMS::String& label) { return self.setPeptideGroupLabel(label); }, "label"_a, "Sets the peptide group label")
        .def("getPeptideGroupLabel", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return self.getPeptideGroupLabel(); }, "Get the peptide group label")
        .def("setChargeState", [](OpenMS::TargetedExperimentHelper::Peptide& self, int charge) { return self.setChargeState(charge); }, "charge"_a, "Sets the peptide or compound charge states")
        .def("hasCharge", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return self.hasCharge(); }, "Whether product has set charge state")
        .def("getChargeState", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return self.getChargeState(); }, "Returns the peptide or compound charge state")
        .def("hasRetentionTime", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return self.hasRetentionTime(); }, "Gets compound or peptide retention time")
        .def("getRetentionTime", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return self.getRetentionTime(); }, "Gets compound or peptide retention time")
        .def("getRetentionTimeType", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return self.getRetentionTimeType(); }, "Get compound or peptide retentiontime type")
        .def("getRetentionTimeUnit", [](const OpenMS::TargetedExperimentHelper::Peptide& self) { return self.getRetentionTimeUnit(); }, "Get compound or peptide retentiontime unit (minute/seconds)")
        .def(nb::self != nb::self)
        
        .def_rw("protein_refs", &OpenMS::TargetedExperimentHelper::Peptide::protein_refs)
        .def_rw("evidence", &OpenMS::TargetedExperimentHelper::Peptide::evidence)
        .def_rw("sequence", &OpenMS::TargetedExperimentHelper::Peptide::sequence)
        .def_rw("mods", &OpenMS::TargetedExperimentHelper::Peptide::mods)
        .def_rw("rts", &OpenMS::TargetedExperimentHelper::Peptide::rts)
        .def_rw("id", &OpenMS::TargetedExperimentHelper::Peptide::id)
        ;
    def_CVTermList<OpenMS::TargetedExperimentHelper::Peptide>(peptide_class);
    def_MetaInfoInterface<OpenMS::TargetedExperimentHelper::Peptide>(peptide_class);

    // -----------------------------------------------------------------------
    // PeptideAndProteinQuant_PeptideData
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideAndProteinQuant::PeptideData>(m, "PeptideAndProteinQuant_PeptideData", "OpenMS class PeptideAndProteinQuant_PeptideData")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::PeptideAndProteinQuant::PeptideData& self) { return OpenMS::PeptideAndProteinQuant::PeptideData(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideAndProteinQuant::PeptideData& self, nb::dict) { return OpenMS::PeptideAndProteinQuant::PeptideData(self); }, "memo"_a)
        .def_rw("abundances", &OpenMS::PeptideAndProteinQuant::PeptideData::abundances)
        .def_rw("psm_counts", &OpenMS::PeptideAndProteinQuant::PeptideData::psm_counts)
        .def_rw("total_abundances", &OpenMS::PeptideAndProteinQuant::PeptideData::total_abundances)
        .def_rw("total_psm_counts", &OpenMS::PeptideAndProteinQuant::PeptideData::total_psm_counts)
        .def_rw("accessions", &OpenMS::PeptideAndProteinQuant::PeptideData::accessions)
        .def_rw("psm_count", &OpenMS::PeptideAndProteinQuant::PeptideData::psm_count)
        ;

    // -----------------------------------------------------------------------
    // PeptideAndProteinQuant_ProteinData
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideAndProteinQuant::ProteinData>(m, "PeptideAndProteinQuant_ProteinData", "OpenMS class PeptideAndProteinQuant_ProteinData")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::PeptideAndProteinQuant::ProteinData& self) { return OpenMS::PeptideAndProteinQuant::ProteinData(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideAndProteinQuant::ProteinData& self, nb::dict) { return OpenMS::PeptideAndProteinQuant::ProteinData(self); }, "memo"_a)
        .def_rw("peptide_abundances", &OpenMS::PeptideAndProteinQuant::ProteinData::peptide_abundances)
        .def_rw("peptide_psm_counts", &OpenMS::PeptideAndProteinQuant::ProteinData::peptide_psm_counts)
        .def_rw("channel_level_abundances", &OpenMS::PeptideAndProteinQuant::ProteinData::channel_level_abundances)
        .def_rw("file_level_psm_counts", &OpenMS::PeptideAndProteinQuant::ProteinData::file_level_psm_counts)
        .def_rw("total_abundances", &OpenMS::PeptideAndProteinQuant::ProteinData::total_abundances)
        .def_rw("total_psm_counts", &OpenMS::PeptideAndProteinQuant::ProteinData::total_psm_counts)
        .def_rw("total_distinct_peptides", &OpenMS::PeptideAndProteinQuant::ProteinData::total_distinct_peptides)
        .def_rw("psm_count", &OpenMS::PeptideAndProteinQuant::ProteinData::psm_count)
        ;

    // -----------------------------------------------------------------------
    // PeptideAndProteinQuant_Statistics
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideAndProteinQuant::Statistics>(m, "PeptideAndProteinQuant_Statistics", "OpenMS class PeptideAndProteinQuant_Statistics")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::PeptideAndProteinQuant::Statistics& self) { return OpenMS::PeptideAndProteinQuant::Statistics(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideAndProteinQuant::Statistics& self, nb::dict) { return OpenMS::PeptideAndProteinQuant::Statistics(self); }, "memo"_a)
        .def_rw("n_samples", &OpenMS::PeptideAndProteinQuant::Statistics::n_samples)
        .def_rw("n_fractions", &OpenMS::PeptideAndProteinQuant::Statistics::n_fractions)
        .def_rw("n_ms_files", &OpenMS::PeptideAndProteinQuant::Statistics::n_ms_files)
        .def_rw("quant_proteins", &OpenMS::PeptideAndProteinQuant::Statistics::quant_proteins)
        .def_rw("too_few_peptides", &OpenMS::PeptideAndProteinQuant::Statistics::too_few_peptides)
        .def_rw("quant_peptides", &OpenMS::PeptideAndProteinQuant::Statistics::quant_peptides)
        .def_rw("total_peptides", &OpenMS::PeptideAndProteinQuant::Statistics::total_peptides)
        .def_rw("quant_features", &OpenMS::PeptideAndProteinQuant::Statistics::quant_features)
        .def_rw("total_features", &OpenMS::PeptideAndProteinQuant::Statistics::total_features)
        .def_rw("blank_features", &OpenMS::PeptideAndProteinQuant::Statistics::blank_features)
        .def_rw("ambig_features", &OpenMS::PeptideAndProteinQuant::Statistics::ambig_features)
        ;

    // -----------------------------------------------------------------------
    // PeptideProteinResolution_ConnectedComponent (ConnectedComponent)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConnectedComponent>(m, "PeptideProteinResolution_ConnectedComponent",
        "Connected component of the bipartite protein-peptide graph")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::ConnectedComponent& self) { return OpenMS::ConnectedComponent(self); })
        .def("__deepcopy__", [](const OpenMS::ConnectedComponent& self, nb::dict) { return OpenMS::ConnectedComponent(self); }, "memo"_a)
        .def_rw("prot_grp_indices", &OpenMS::ConnectedComponent::prot_grp_indices)
        .def_rw("pep_indices", &OpenMS::ConnectedComponent::pep_indices)
        ;

    // -----------------------------------------------------------------------
    // PeptideProteinResolution
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideProteinResolution>(m, "PeptideProteinResolution", 
        R"doc(
Resolves shared peptides based on protein scores
Resolves connected components of the bipartite protein-peptide graph based
on protein probabilities/scores and adds them as additional protein_groups
to the protein identification run processed.
Thereby greedily assigns shared peptides in this component uniquely to the
proteins of the current @em best @em indistinguishable protein group, until
every peptide is uniquely assigned. This effectively allows more peptides to
be used in ProteinQuantifier at the cost of potentially additional noise in
the peptides quantities.
In accordance with most state-of-the-art protein inference tools, only the
best hit (PSM) for a peptide ID is considered.  Probability ties are
currently resolved by taking the protein with larger number of peptides
The class could provide iterator for ConnectedComponents in the
future. One could extend the graph to include all PeptideHits (not only the
best). It becomes a tripartite graph with larger connected components then.
Maybe extend it to work with MS1 features. Separate resolution and adding
groups to output
)doc")
        .def(nb::init<bool>())
        .def("__copy__", [](const OpenMS::PeptideProteinResolution& self) { return OpenMS::PeptideProteinResolution(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideProteinResolution& self, nb::dict) { return OpenMS::PeptideProteinResolution(self); }, "memo"_a)
        .def_static("run", [](std::vector<OpenMS::ProteinIdentification> inferred_protein_id, OpenMS::PeptideIdentificationList& inferred_peptide_ids) {
            OpenMS::PeptideProteinResolution::run(inferred_protein_id, inferred_peptide_ids);
            return inferred_protein_id;
        }, "inferred_protein_id"_a, "inferred_peptide_ids"_a,
            R"doc(
Resolves connected components based on posterior probabilities and adds them
as additional protein_groups to the output idXML.
Thereby greedily assigns shared peptides in this component uniquely to
the proteins of the current BEST INDISTINGUISHABLE protein group,
ready to be used in ProteinQuantifier then.
This is achieved by removing all other evidence from the input
PeptideIDs and iterating until each peptide is uniquely assigned.
In accordance with Fido only the best hit (PSM) for an ID is considered.
Probability ties resolved by taking protein with largest number of peptides
:param conn_comp: The component to be resolved
:param protein: ProteinIdentification object storing IDs and groups
:param peptides: Vector of ProteinIdentifications with links to the proteins
:returns: Updated inferred_protein_id
static members
)doc")
        .def("buildGraph", [](OpenMS::PeptideProteinResolution& self, OpenMS::ProteinIdentification& protein, const OpenMS::PeptideIdentificationList& peptides, bool skip_sort) { return self.buildGraph(protein, peptides, skip_sort); }, "protein"_a, "peptides"_a, "skip_sort"_a = false)
        .def("resolveGraph", [](OpenMS::PeptideProteinResolution& self, OpenMS::ProteinIdentification& protein, OpenMS::PeptideIdentificationList& peptides) { return self.resolveGraph(protein, peptides); }, "protein"_a, "peptides"_a, 
            R"doc(
Initialize and store the graph (= maps), needs sorted groups for
correct functionality. Therefore sorts the indist. protein groups
if not skipped
:param protein: ProteinIdentification object storing IDs and groups
:param peptides: Vector of ProteinIdentifications with links to the proteins
:param skip_sort: Skips sorting of groups, nothing is modified then
)doc")
        .def("findConnectedComponent", [](OpenMS::PeptideProteinResolution& self, size_t root_prot_grp) { return self.findConnectedComponent(root_prot_grp); }, "root_prot_grp"_a,
            R"doc(
Applies resolveConnectedComponent to every component of the graph and
is able to write statistics when specified. Parameters will
both be mutated in this method
:param protein: ProteinIdentification object storing IDs and groups
:param peptides: vector of ProteinIdentifications with links to the proteins
)doc")
        .def("resolveConnectedComponent", [](OpenMS::PeptideProteinResolution& self, OpenMS::ConnectedComponent& conn_comp, OpenMS::ProteinIdentification& protein, OpenMS::PeptideIdentificationList& peptides) { return self.resolveConnectedComponent(conn_comp, protein, peptides); }, "conn_comp"_a, "protein"_a, "peptides"_a, 
            R"doc(
Does a BFS on the two maps (= two parts of the graph; indist. prot. groups
and peptides), switching from one to the other in each step
:param root_prot_grp: Starts the BFS at this protein group index
:return: Returns a Connected Component as set of group and peptide indices
)doc")
        ;

    // -----------------------------------------------------------------------
    // PercolatorFeatureSetHelper
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PercolatorFeatureSetHelper>(m, "PercolatorFeatureSetHelper", 
        R"doc(
Percolator feature set and integration helper
This class contains functions to handle (compute, aggregate, integrate)
Percolator features. This includes the calculation or extraction of
Percolator features depending on the search engine(s) for later use with
PercolatorAdapter. It also includes handling the reintegration of the
percolator result into the set of Identifications
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PercolatorFeatureSetHelper &>())
        .def("__copy__", [](const OpenMS::PercolatorFeatureSetHelper& self) { return OpenMS::PercolatorFeatureSetHelper(self); })
        .def("__deepcopy__", [](const OpenMS::PercolatorFeatureSetHelper& self, nb::dict) { return OpenMS::PercolatorFeatureSetHelper(self); }, "memo"_a)
        .def_static("concatMULTISEPeptideIds", [](OpenMS::PeptideIdentificationList& all_peptide_ids, OpenMS::PeptideIdentificationList& new_peptide_ids, const OpenMS::String& search_engine) { return OpenMS::PercolatorFeatureSetHelper::concatMULTISEPeptideIds(all_peptide_ids, new_peptide_ids, search_engine); }, "all_peptide_ids"_a, "new_peptide_ids"_a, "search_engine"_a)
        .def_static("mergeMULTISEPeptideIds", [](OpenMS::PeptideIdentificationList& all_peptide_ids, OpenMS::PeptideIdentificationList& new_peptide_ids, const OpenMS::String& search_engine) { return OpenMS::PercolatorFeatureSetHelper::mergeMULTISEPeptideIds(all_peptide_ids, new_peptide_ids, search_engine); }, "all_peptide_ids"_a, "new_peptide_ids"_a, "search_engine"_a, 
            R"doc(
Appends a vector of PeptideIdentification to another and prepares Percolator features in MetaInfo (With the respective key "CONCAT:" + search_engine)
:param all_peptide_ids: PeptideIdentification vector to append to
:param new_peptide_ids: PeptideIdentification vector to be appended
:param search_engine: Search engine to depend on for feature creation
)doc")
        .def_static("mergeMULTISEProteinIds", [](std::vector<OpenMS::ProteinIdentification> all_protein_ids, std::vector<OpenMS::ProteinIdentification> new_protein_ids) {
            OpenMS::PercolatorFeatureSetHelper::mergeMULTISEProteinIds(all_protein_ids, new_protein_ids);
            return nb::make_tuple(all_protein_ids, new_protein_ids);
        }, "all_protein_ids"_a, "new_protein_ids"_a,
            R"doc(
Merges a vector of PeptideIdentification into another and prepares the merged MetaInfo and scores for collection in addMULTISEFeatures for feature registration
:param all_peptide_idsL: PeptideIdentification vector to be merged into
:param new_peptide_idsL: PeptideIdentification vector to merge
:param search_engineL: Search engine to create features from their scores
:returns: Tuple of (updated all_protein_ids, updated new_protein_ids)
)doc")
        .def_static("addMSGFFeatures", [](OpenMS::PeptideIdentificationList& peptide_ids, std::vector<OpenMS::String> feature_set) {
            OpenMS::PercolatorFeatureSetHelper::addMSGFFeatures(peptide_ids, feature_set);
            return feature_set;
        }, "peptide_ids"_a, "feature_set"_a,
            R"doc(
Concatenates SearchParameter of multiple search engine runs and merges PeptideEvidences, collects used search engines in MetaInfo for collection in addMULTISEFeatures for feature registration
:param all_protein_ids: ProteinIdentification vector to be merged into
:param new_protein_ids: ProteinIdentification vector to merge
:returns: Updated feature_set
)doc")
        .def_static("addXTANDEMFeatures", [](OpenMS::PeptideIdentificationList& peptide_ids, std::vector<OpenMS::String> feature_set) {
            OpenMS::PercolatorFeatureSetHelper::addXTANDEMFeatures(peptide_ids, feature_set);
            return feature_set;
        }, "peptide_ids"_a, "feature_set"_a,
            R"doc(
Creates and adds MSGF+ specific Percolator features and registers them in feature_set. MSGF+ should be run with the addFeatures flag enabled
:param peptide_ids: PeptideIdentification vector to create Percolator features in
:param feature_set: Register of added features
:returns: Updated feature_set
)doc")
        .def_static("addCOMETFeatures", [](OpenMS::PeptideIdentificationList& peptide_ids, std::vector<OpenMS::String> feature_set) {
            OpenMS::PercolatorFeatureSetHelper::addCOMETFeatures(peptide_ids, feature_set);
            return feature_set;
        }, "peptide_ids"_a, "feature_set"_a,
            R"doc(
Creates and adds X!Tandem specific Percolator features and registers them in feature_set
:param peptide_ids: PeptideIdentification vector to create Percolator features in
:param feature_set: Register of added features
:returns: Updated feature_set
)doc")
        .def_static("addMASCOTFeatures", [](OpenMS::PeptideIdentificationList& peptide_ids, std::vector<OpenMS::String> feature_set) {
            OpenMS::PercolatorFeatureSetHelper::addMASCOTFeatures(peptide_ids, feature_set);
            return feature_set;
        }, "peptide_ids"_a, "feature_set"_a,
            R"doc(
Creates and adds Comet specific Percolator features and registers them in feature_set
:param peptide_ids: PeptideIdentification vector to create Percolator features in
:param feature_set: Register of added features
:returns: Updated feature_set
)doc")
        .def_static("addMULTISEFeatures", [](OpenMS::PeptideIdentificationList& peptide_ids, std::vector<OpenMS::String> search_engines_used, std::vector<OpenMS::String> feature_set, bool complete_only, bool limits_imputation) {
            OpenMS::PercolatorFeatureSetHelper::addMULTISEFeatures(peptide_ids, search_engines_used, feature_set, complete_only, limits_imputation);
            return feature_set;
        }, "peptide_ids"_a, "search_engines_used"_a, "feature_set"_a, "complete_only"_a, "limits_imputation"_a,
            R"doc(
Creates and adds Mascot specific Percolator features and registers them in feature_set
:param peptide_ids: PeptideIdentification vector to create Percolator features in
:param feature_set: Register of added features
:returns: Updated feature_set
)doc")
        .def_static("addCONCATSEFeatures", [](OpenMS::PeptideIdentificationList& peptide_id_list, std::vector<OpenMS::String> search_engines_used, std::vector<OpenMS::String> feature_set) {
            OpenMS::PercolatorFeatureSetHelper::addCONCATSEFeatures(peptide_id_list, search_engines_used, feature_set);
            return feature_set;
        }, "peptide_id_list"_a, "search_engines_used"_a, "feature_set"_a,
            R"doc(
Adds multiple search engine specific Percolator features and registers them in feature_set
:param peptide_ids: PeptideIdentification vector to create Percolator features in
:param search_engines_used: The list of search engines to be considered
:param feature_set: Register of added features
:param complete_only: Will only add features for PeptideIdentifications where all given search engines identified something
:param limits_imputation: Uses C++ numeric limits as imputed values instead of min/max of that feature
:returns: Updated feature_set
)doc")
        .def_static("checkExtraFeatures", [](const std::vector<OpenMS::PeptideHit>& psms, std::vector<OpenMS::String> extra_features) {
            OpenMS::PercolatorFeatureSetHelper::checkExtraFeatures(psms, extra_features);
            return extra_features;
        }, "psms"_a, "extra_features"_a,
            R"doc(
Adds multiple search engine specific Percolator features and registers them in feature_set
This struct can be used to store both peak or feature indices
:param peptide_ids: PeptideIdentification vector to create Percolator features in
:param search_engines_used: The list of search engines to be considered
:param feature_set: Register of added features
:returns: Updated extra_features with unavailable features removed
)doc")
        ;

    // -----------------------------------------------------------------------
    // PrecalAveragine
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHHelperClasses::PrecalculatedAveragine>(m, "PrecalAveragine", 
        R"doc(
Averagine patterns pre-calculated for speed up.
Used for fast isotope cosine calculation.
Constructors
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHHelperClasses::PrecalculatedAveragine &>())
        .def("__copy__", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self) { return OpenMS::FLASHHelperClasses::PrecalculatedAveragine(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, nb::dict) { return OpenMS::FLASHHelperClasses::PrecalculatedAveragine(self); }, "memo"_a)
        .def("get", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.get(mass); }, "mass"_a, "Get isotope distribution for given mass")
        .def("getMaxIsotopeIndex", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self) { return self.getMaxIsotopeIndex(); }, "Get max isotope index")
        .def("setMaxIsotopeIndex", [](OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, int index) { return self.setMaxIsotopeIndex(index); }, "index"_a, "Set max isotope index")
        .def("getLeftCountFromApex", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.getLeftCountFromApex(mass); }, "mass"_a, "Get isotope count left of apex")
        .def("getRightCountFromApex", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.getRightCountFromApex(mass); }, "mass"_a, "Get isotope count right of apex")
        .def("getApexIndex", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.getApexIndex(mass); }, "mass"_a, "Get apex isotope index")
        .def("getLastIndex", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.getLastIndex(mass); }, "mass"_a, "Get last isotope index")
        .def("getAverageMassDelta", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.getAverageMassDelta(mass); }, "mass"_a, "Get mass diff between avg and mono")
        .def("getMostAbundantMassDelta", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.getMostAbundantMassDelta(mass); }, "mass"_a, "Get mass diff between most abundant and mono")
        .def("getSNRMultiplicationFactor", [](const OpenMS::FLASHHelperClasses::PrecalculatedAveragine& self, double mass) { return self.getSNRMultiplicationFactor(mass); }, "mass"_a, "Get SNR multiplication factor")

        .def("__init__", [](OpenMS::FLASHHelperClasses::PrecalculatedAveragine* self, double min_mass, double max_mass, double delta, OpenMS::CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine) {
            new (self) OpenMS::FLASHHelperClasses::PrecalculatedAveragine(min_mass, max_mass, delta, generator, use_RNA_averagine);
        }, "min_mass"_a, "max_mass"_a, "delta"_a, "generator"_a, "use_RNA_averagine"_a)

        .def("__init__", [](OpenMS::FLASHHelperClasses::PrecalculatedAveragine* self, double min_mass, double max_mass, double delta, OpenMS::CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine, double decoy_iso_distance) {
            new (self) OpenMS::FLASHHelperClasses::PrecalculatedAveragine(min_mass, max_mass, delta, generator, use_RNA_averagine, decoy_iso_distance);
        }, "min_mass"_a, "max_mass"_a, "delta"_a, "generator"_a, "use_RNA_averagine"_a, "decoy_iso_distance"_a)
        ;

    // -----------------------------------------------------------------------
    // PrecursorPurity
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PrecursorPurity>(m, "PrecursorPurity", 
        R"doc(
Precursor purity or noise estimation
This class computes metrics for precursor isolation window purity (or noise)
The function extracts the peaks from an isolation window targeted for fragmentation
and determines which peaks are isotopes of the target and which come from other sources
The intensities of the assumed target peaks are summed up as the target intensity
Using this information it calculates an intensity ratio for the relative intensity of the target
compared to other sources
These metrics are combined over the previous and the next MS1 spectrum
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PrecursorPurity &>())
        .def("__copy__", [](const OpenMS::PrecursorPurity& self) { return OpenMS::PrecursorPurity(self); })
        .def("__deepcopy__", [](const OpenMS::PrecursorPurity& self, nb::dict) { return OpenMS::PrecursorPurity(self); }, "memo"_a)
        .def_static("computePrecursorPurity", [](const OpenMS::MSSpectrum& ms1, const OpenMS::Precursor& pre, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm) { return OpenMS::PrecursorPurity::computePrecursorPurity(ms1, pre, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm); }, "ms1"_a, "pre"_a, "precursor_mass_tolerance"_a, "precursor_mass_tolerance_unit_ppm"_a, 
            R"doc(
)doc")
        ;

    // -----------------------------------------------------------------------
    // PurityScores (PrecursorPurity::PurityScores)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PrecursorPurity::PurityScores>(m, "PurityScores",
        "Precursor purity scores for a single precursor")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::PrecursorPurity::PurityScores& self) { return OpenMS::PrecursorPurity::PurityScores(self); })
        .def("__deepcopy__", [](const OpenMS::PrecursorPurity::PurityScores& self, nb::dict) { return OpenMS::PrecursorPurity::PurityScores(self); }, "memo"_a)
        .def_rw("total_intensity", &OpenMS::PrecursorPurity::PurityScores::total_intensity)
        .def_rw("target_intensity", &OpenMS::PrecursorPurity::PurityScores::target_intensity)
        .def_rw("signal_proportion", &OpenMS::PrecursorPurity::PurityScores::signal_proportion)
        .def_rw("target_peak_count", &OpenMS::PrecursorPurity::PurityScores::target_peak_count)
        .def_rw("interfering_peak_count", &OpenMS::PrecursorPurity::PurityScores::interfering_peak_count)
        .def_rw("interfering_peaks", &OpenMS::PrecursorPurity::PurityScores::interfering_peaks)
        ;

    // -----------------------------------------------------------------------
    // ProbablePhosphoSites
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProbablePhosphoSites>(m, "ProbablePhosphoSites", "DefaultParamHandler")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProbablePhosphoSites &>())
        .def("__copy__", [](const OpenMS::ProbablePhosphoSites& self) { return OpenMS::ProbablePhosphoSites(self); })
        .def("__deepcopy__", [](const OpenMS::ProbablePhosphoSites& self, nb::dict) { return OpenMS::ProbablePhosphoSites(self); }, "memo"_a)
        .def_rw("first", &OpenMS::ProbablePhosphoSites::first)
        .def_rw("second", &OpenMS::ProbablePhosphoSites::second)
        .def_rw("seq_1", &OpenMS::ProbablePhosphoSites::seq_1)
        .def_rw("seq_2", &OpenMS::ProbablePhosphoSites::seq_2)
        .def_rw("peak_depth", &OpenMS::ProbablePhosphoSites::peak_depth)
        .def_rw("AScore", &OpenMS::ProbablePhosphoSites::AScore)
        ;

    // -----------------------------------------------------------------------
    // ProteinInference
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteinInference>(m, "ProteinInference", 
        R"doc(
[experimental class] given a peptide quantitation, infer corresponding protein quantities
Infers protein ratios from peptide ratios (currently using unique peptides only).
Use the IDMapper class to add protein and peptide information to a
quantitative ConsensusMap prior to this step
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProteinInference &>())
        .def("__copy__", [](const OpenMS::ProteinInference& self) { return OpenMS::ProteinInference(self); })
        .def("__deepcopy__", [](const OpenMS::ProteinInference& self, nb::dict) { return OpenMS::ProteinInference(self); }, "memo"_a)
        .def("infer", [](OpenMS::ProteinInference& self, OpenMS::ConsensusMap& consensus_map, unsigned int reference_map) { return self.infer(consensus_map, reference_map); }, "consensus_map"_a, "reference_map"_a)
        ;

    // -----------------------------------------------------------------------
    // Scores
    // -----------------------------------------------------------------------
    auto scores_class = nb::class_<OpenMS::Scores>(m, "Scores", 
        R"doc(
Utility class for score type handling in identification and quantification workflows.
This class provides centralized handling of score types used in peptide/protein
identification, quantification, and PTM localization. It defines the hierarchy of
score types and provides utility methods for score type conversion, comparison, and lookup.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Scores &>())
        .def("__copy__", [](const OpenMS::Scores& self) { return OpenMS::Scores(self); })
        .def("__deepcopy__", [](const OpenMS::Scores& self, nb::dict) { return OpenMS::Scores(self); }, "memo"_a)

        .def_static("isScoreType", &OpenMS::Scores::isScoreType, "score_name"_a, "type"_a,
            "Check if the given score name corresponds to a specific ID score type")

        .def_static("parseIDType", &OpenMS::Scores::parseIDType, "score_type"_a,
            "Convert a string representation of an ID score type to an IDType enum")

        .def_static("isHigherBetter", &OpenMS::Scores::isHigherBetter, "type"_a,
            "Determine whether a higher score is better for the given ID score type")

        .def_static("getAllIDScoreNames", &OpenMS::Scores::getAllIDScoreNames,
            "Get a vector of all ID score names used in OpenMS")

        .def_static("normalizeScoreName", &OpenMS::Scores::normalizeScoreName, "score_name"_a,
            "Normalize a score name by removing the '_score' suffix if present")

        .def_static("isKnownScoreType", &OpenMS::Scores::isKnownScoreType, "score_name"_a,
            "Check if a score name is a known score type after normalization")
        ;
    // IDType enum nested under Scores
    nb::enum_<OpenMS::Scores::IDType>(scores_class, "IDType", "Enum for identification score types", nb::is_arithmetic())
        .value("RAW", OpenMS::Scores::IDType::RAW)
        .value("RAW_EVAL", OpenMS::Scores::IDType::RAW_EVAL)
        .value("PP", OpenMS::Scores::IDType::PP)
        .value("PEP", OpenMS::Scores::IDType::PEP)
        .value("FDR", OpenMS::Scores::IDType::FDR)
        .value("QVAL", OpenMS::Scores::IDType::QVAL)
        ;

    // -----------------------------------------------------------------------
    // SelectorParameters
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeatureSelector::SelectorParameters>(m, "SelectorParameters", 
        R"doc(
Structure to configure MRMFeatureSelector parameters
Attributes:
- nn_threshold: Number of nearest neighbors to include in optimization (default: 4)
- locality_weight: Weight compounds by retention time proximity (default: False)
- select_transition_group: Use component groups instead of components (default: True)
- segment_window_length: Number of components in the sliding window (default: 8)
- segment_step_length: Step size for sliding window (default: 4)
- optimal_threshold: Cutoff for selection, 0 < x < 1 (default: 0.5)
Example:
.. code-block:: python
from pyopenms import *
params = SelectorParameters()
params.nn_threshold = 6
params.optimal_threshold = 0.7
params.select_transition_group = False
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeatureSelector::SelectorParameters &>())
        .def("__copy__", [](const OpenMS::MRMFeatureSelector::SelectorParameters& self) { return OpenMS::MRMFeatureSelector::SelectorParameters(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeatureSelector::SelectorParameters& self, nb::dict) { return OpenMS::MRMFeatureSelector::SelectorParameters(self); }, "memo"_a)
        .def_rw("nn_threshold", &OpenMS::MRMFeatureSelector::SelectorParameters::nn_threshold)
        .def_rw("locality_weight", &OpenMS::MRMFeatureSelector::SelectorParameters::locality_weight)
        .def_rw("select_transition_group", &OpenMS::MRMFeatureSelector::SelectorParameters::select_transition_group)
        .def_rw("segment_window_length", &OpenMS::MRMFeatureSelector::SelectorParameters::segment_window_length)
        .def_rw("segment_step_length", &OpenMS::MRMFeatureSelector::SelectorParameters::segment_step_length)
        .def_rw("optimal_threshold", &OpenMS::MRMFeatureSelector::SelectorParameters::optimal_threshold)
        ;

    // -----------------------------------------------------------------------
    // SimpleOpenMSSpectraFactory
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SimpleOpenMSSpectraFactory>(m, "SimpleOpenMSSpectraFactory", "A factory method that returns two ISpectrumAccess implementations")
        .def(nb::init<>())
        ;

    // -----------------------------------------------------------------------
    // SiriusMSFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SiriusMSFile>(m, "SiriusMSFile", "OpenMS class SiriusMSFile")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SiriusMSFile &>())
        .def("__copy__", [](const OpenMS::SiriusMSFile& self) { return OpenMS::SiriusMSFile(self); })
        .def("__deepcopy__", [](const OpenMS::SiriusMSFile& self, nb::dict) { return OpenMS::SiriusMSFile(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // SiriusMSFile_AccessionInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SiriusMSFile::AccessionInfo>(m, "SiriusMSFile_AccessionInfo", "OpenMS class SiriusMSFile_AccessionInfo")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SiriusMSFile::AccessionInfo &>())
        .def("__copy__", [](const OpenMS::SiriusMSFile::AccessionInfo& self) { return OpenMS::SiriusMSFile::AccessionInfo(self); })
        .def("__deepcopy__", [](const OpenMS::SiriusMSFile::AccessionInfo& self, nb::dict) { return OpenMS::SiriusMSFile::AccessionInfo(self); }, "memo"_a)
        .def_rw("sf_path", &OpenMS::SiriusMSFile::AccessionInfo::sf_path)
        .def_rw("sf_type", &OpenMS::SiriusMSFile::AccessionInfo::sf_type)
        .def_rw("sf_filename", &OpenMS::SiriusMSFile::AccessionInfo::sf_filename)
        .def_rw("sf_accession", &OpenMS::SiriusMSFile::AccessionInfo::sf_accession)
        .def_rw("native_id_accession", &OpenMS::SiriusMSFile::AccessionInfo::native_id_accession)
        .def_rw("native_id_type", &OpenMS::SiriusMSFile::AccessionInfo::native_id_type)
        ;

    // -----------------------------------------------------------------------
    // SiriusMSFile_CompoundInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SiriusMSFile::CompoundInfo>(m, "SiriusMSFile_CompoundInfo", "OpenMS class SiriusMSFile_CompoundInfo")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SiriusMSFile::CompoundInfo &>())
        .def("__copy__", [](const OpenMS::SiriusMSFile::CompoundInfo& self) { return OpenMS::SiriusMSFile::CompoundInfo(self); })
        .def("__deepcopy__", [](const OpenMS::SiriusMSFile::CompoundInfo& self, nb::dict) { return OpenMS::SiriusMSFile::CompoundInfo(self); }, "memo"_a)
        .def_rw("cmp", &OpenMS::SiriusMSFile::CompoundInfo::cmp)
        .def_rw("pmass", &OpenMS::SiriusMSFile::CompoundInfo::pmass)
        .def_rw("pint_mono", &OpenMS::SiriusMSFile::CompoundInfo::pint_mono)
        .def_rw("rt", &OpenMS::SiriusMSFile::CompoundInfo::rt)
        .def_rw("fmz", &OpenMS::SiriusMSFile::CompoundInfo::fmz)
        .def_rw("fid", &OpenMS::SiriusMSFile::CompoundInfo::fid)
        .def_rw("formula", &OpenMS::SiriusMSFile::CompoundInfo::formula)
        .def_rw("charge", &OpenMS::SiriusMSFile::CompoundInfo::charge)
        .def_rw("ionization", &OpenMS::SiriusMSFile::CompoundInfo::ionization)
        .def_rw("des", &OpenMS::SiriusMSFile::CompoundInfo::des)
        .def_rw("specref_format", &OpenMS::SiriusMSFile::CompoundInfo::specref_format)
        .def_rw("source_file", &OpenMS::SiriusMSFile::CompoundInfo::source_file)
        .def_rw("source_format", &OpenMS::SiriusMSFile::CompoundInfo::source_format)
        .def_rw("native_ids", &OpenMS::SiriusMSFile::CompoundInfo::native_ids)
        .def_rw("native_ids_id", &OpenMS::SiriusMSFile::CompoundInfo::native_ids_id)
        .def_rw("m_ids", &OpenMS::SiriusMSFile::CompoundInfo::m_ids)
        .def_rw("m_ids_id", &OpenMS::SiriusMSFile::CompoundInfo::m_ids_id)
        .def_rw("scan_indices", &OpenMS::SiriusMSFile::CompoundInfo::scan_indices)
        .def_rw("specrefs", &OpenMS::SiriusMSFile::CompoundInfo::specrefs)
        .def_rw("file_index", &OpenMS::SiriusMSFile::CompoundInfo::file_index)
        ;

    // -----------------------------------------------------------------------
    // SpectralMatch
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectralMatch>(m, "SpectralMatch", "OpenMS class SpectralMatch")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectralMatch &>())
        .def("__copy__", [](const OpenMS::SpectralMatch& self) { return OpenMS::SpectralMatch(self); })
        .def("__deepcopy__", [](const OpenMS::SpectralMatch& self, nb::dict) { return OpenMS::SpectralMatch(self); }, "memo"_a)
        .def("getObservedPrecursorMass", [](const OpenMS::SpectralMatch& self) { return self.getObservedPrecursorMass(); })
        .def("setObservedPrecursorMass", [](OpenMS::SpectralMatch& self, const double& p0) { return self.setObservedPrecursorMass(p0); })
        .def("getObservedPrecursorRT", [](const OpenMS::SpectralMatch& self) { return self.getObservedPrecursorRT(); })
        .def("setObservedPrecursorRT", [](OpenMS::SpectralMatch& self, const double& p0) { return self.setObservedPrecursorRT(p0); })
        .def("getFoundPrecursorMass", [](const OpenMS::SpectralMatch& self) { return self.getFoundPrecursorMass(); })
        .def("setFoundPrecursorMass", [](OpenMS::SpectralMatch& self, const double& p0) { return self.setFoundPrecursorMass(p0); })
        .def("getFoundPrecursorCharge", [](const OpenMS::SpectralMatch& self) { return self.getFoundPrecursorCharge(); })
        .def("setFoundPrecursorCharge", [](OpenMS::SpectralMatch& self, const int& p0) { return self.setFoundPrecursorCharge(p0); })
        .def("getMatchingScore", [](const OpenMS::SpectralMatch& self) { return self.getMatchingScore(); })
        .def("setMatchingScore", [](OpenMS::SpectralMatch& self, const double& p0) { return self.setMatchingScore(p0); })
        .def("getObservedSpectrumIndex", [](const OpenMS::SpectralMatch& self) { return self.getObservedSpectrumIndex(); })
        .def("setObservedSpectrumIndex", [](OpenMS::SpectralMatch& self, const size_t& p0) { return self.setObservedSpectrumIndex(p0); })
        .def("getMatchingSpectrumIndex", [](const OpenMS::SpectralMatch& self) { return self.getMatchingSpectrumIndex(); })
        .def("setMatchingSpectrumIndex", [](OpenMS::SpectralMatch& self, const size_t& p0) { return self.setMatchingSpectrumIndex(p0); })
        .def("getPrimaryIdentifier", [](const OpenMS::SpectralMatch& self) { return self.getPrimaryIdentifier(); })
        .def("setPrimaryIdentifier", [](OpenMS::SpectralMatch& self, const OpenMS::String& p0) { return self.setPrimaryIdentifier(p0); })
        .def("getSecondaryIdentifier", [](const OpenMS::SpectralMatch& self) { return self.getSecondaryIdentifier(); })
        .def("setSecondaryIdentifier", [](OpenMS::SpectralMatch& self, const OpenMS::String& p0) { return self.setSecondaryIdentifier(p0); })
        .def("getCommonName", [](const OpenMS::SpectralMatch& self) { return self.getCommonName(); })
        .def("setCommonName", [](OpenMS::SpectralMatch& self, const OpenMS::String& p0) { return self.setCommonName(p0); })
        .def("getSumFormula", [](const OpenMS::SpectralMatch& self) { return self.getSumFormula(); })
        .def("setSumFormula", [](OpenMS::SpectralMatch& self, const OpenMS::String& p0) { return self.setSumFormula(p0); })
        .def("getInchiString", [](const OpenMS::SpectralMatch& self) { return self.getInchiString(); })
        .def("setInchiString", [](OpenMS::SpectralMatch& self, const OpenMS::String& p0) { return self.setInchiString(p0); })
        .def("getSMILESString", [](const OpenMS::SpectralMatch& self) { return self.getSMILESString(); })
        .def("setSMILESString", [](OpenMS::SpectralMatch& self, const OpenMS::String& p0) { return self.setSMILESString(p0); })
        .def("getPrecursorAdduct", [](const OpenMS::SpectralMatch& self) { return self.getPrecursorAdduct(); })
        .def("setPrecursorAdduct", [](OpenMS::SpectralMatch& self, const OpenMS::String& p0) { return self.setPrecursorAdduct(p0); })
        ;

    // -----------------------------------------------------------------------
    // SwathWindowLoader
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SwathWindowLoader>(m, "SwathWindowLoader", 
        R"doc(
Class to read a file describing the Swath Windows * * The file must
of be tab delimited and of the following format: * window_lower
window_upper * 400 425 * 425 450 * ... * * Note that the first line is
a header and will be skipped. *
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SwathWindowLoader &>())
        .def("__copy__", [](const OpenMS::SwathWindowLoader& self) { return OpenMS::SwathWindowLoader(self); })
        .def("__deepcopy__", [](const OpenMS::SwathWindowLoader& self, nb::dict) { return OpenMS::SwathWindowLoader(self); }, "memo"_a)
        .def_static("annotateSwathMapsFromFile", [](const OpenMS::String& filename, std::vector<OpenSwath::SwathMap> swath_maps, bool do_sort, bool force) {
            OpenMS::SwathWindowLoader::annotateSwathMapsFromFile(filename, swath_maps, do_sort, force);
            return swath_maps;
        }, "filename"_a, "swath_maps"_a, "do_sort"_a, "force"_a)
        .def_static("readSwathWindows", [](const OpenMS::String& filename) { std::vector<double> swath_prec_lower; std::vector<double> swath_prec_upper; OpenMS::SwathWindowLoader::readSwathWindows(filename, swath_prec_lower, swath_prec_upper); return std::make_tuple(swath_prec_lower, swath_prec_upper); }, "filename"_a)
        ;

    // -----------------------------------------------------------------------
    // TargetedExperiment
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperiment>(m, "TargetedExperiment", 
        R"doc(
A description of a targeted experiment containing precursor and
production ions
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TargetedExperiment &>())
        .def("__copy__", [](const OpenMS::TargetedExperiment& self) { return OpenMS::TargetedExperiment(self); })
        .def("__deepcopy__", [](const OpenMS::TargetedExperiment& self, nb::dict) { return OpenMS::TargetedExperiment(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self + nb::self)
        .def("clear", [](OpenMS::TargetedExperiment& self, bool clear_meta_data) { return self.clear(clear_meta_data); }, "clear_meta_data"_a)
        .def("setCVs", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::TargetedExperimentHelper::CV>& cvs) { return self.setCVs(cvs); }, "cvs"_a)
        .def("getCVs", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::TargetedExperimentHelper::CV> & { return self.getCVs(); }, nb::rv_policy::reference_internal)
        .def("addCV", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperimentHelper::CV& cv) { return self.addCV(cv); }, "cv"_a)
        .def("setTargetCVTerms", [](OpenMS::TargetedExperiment& self, const OpenMS::CVTermList& cv_terms) { return self.setTargetCVTerms(cv_terms); }, "cv_terms"_a)
        .def("getTargetCVTerms", [](const OpenMS::TargetedExperiment& self) -> const OpenMS::CVTermList & { return self.getTargetCVTerms(); }, nb::rv_policy::reference_internal)
        .def("addTargetCVTerm", [](OpenMS::TargetedExperiment& self, const OpenMS::CVTerm& cv_term) { return self.addTargetCVTerm(cv_term); }, "cv_term"_a)
        .def("setTargetMetaValue", [](OpenMS::TargetedExperiment& self, const OpenMS::String& name, const OpenMS::DataValue& value) { return self.setTargetMetaValue(name, value); }, "name"_a, "value"_a)
        .def("setContacts", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::TargetedExperimentHelper::Contact>& contacts) { return self.setContacts(contacts); }, "contacts"_a)
        .def("getContacts", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Contact> & { return self.getContacts(); }, nb::rv_policy::reference_internal)
        .def("addContact", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperimentHelper::Contact& contact) { return self.addContact(contact); }, "contact"_a)
        .def("setPublications", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::TargetedExperimentHelper::Publication>& publications) { return self.setPublications(publications); }, "publications"_a)
        .def("getPublications", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Publication> & { return self.getPublications(); }, nb::rv_policy::reference_internal)
        .def("addPublication", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperimentHelper::Publication& publication) { return self.addPublication(publication); }, "publication"_a)
        .def("setInstruments", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::TargetedExperimentHelper::Instrument>& instruments) { return self.setInstruments(instruments); }, "instruments"_a)
        .def("getInstruments", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Instrument> & { return self.getInstruments(); }, nb::rv_policy::reference_internal)
        .def("addInstrument", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperimentHelper::Instrument& instrument) { return self.addInstrument(instrument); }, "instrument"_a)
        .def("setSoftware", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::Software>& software) { return self.setSoftware(software); }, "software"_a)
        .def("getSoftware", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::Software> & { return self.getSoftware(); }, nb::rv_policy::reference_internal)
        .def("addSoftware", [](OpenMS::TargetedExperiment& self, const OpenMS::Software& software) { return self.addSoftware(software); }, "software"_a)
        .def("setProteins", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::TargetedExperimentHelper::Protein>& proteins) { return self.setProteins(proteins); }, "proteins"_a)
        .def("setProteins", [](OpenMS::TargetedExperiment& self, std::vector<OpenMS::TargetedExperimentHelper::Protein>& proteins) { return self.setProteins(proteins); }, "proteins"_a)
        .def("getProteins", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Protein> & { return self.getProteins(); }, nb::rv_policy::reference_internal)
        .def("getProteinByRef", [](const OpenMS::TargetedExperiment& self, const OpenMS::String& ref) -> const OpenMS::TargetedExperimentHelper::Protein & { return self.getProteinByRef(ref); }, "ref"_a, nb::rv_policy::reference_internal)
        .def("hasProtein", [](const OpenMS::TargetedExperiment& self, const OpenMS::String& ref) { return self.hasProtein(ref); }, "ref"_a)
        .def("addProtein", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperimentHelper::Protein& protein) { return self.addProtein(protein); }, "protein"_a)
        .def("setCompounds", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::TargetedExperimentHelper::Compound>& rhs) { return self.setCompounds(rhs); }, "rhs"_a)
        .def("getCompounds", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Compound> & { return self.getCompounds(); }, nb::rv_policy::reference_internal)
        .def("addCompound", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperimentHelper::Compound& rhs) { return self.addCompound(rhs); }, "rhs"_a)
        .def("setPeptides", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::TargetedExperimentHelper::Peptide>& rhs) { return self.setPeptides(rhs); }, "rhs"_a)
        .def("setPeptides", [](OpenMS::TargetedExperiment& self, std::vector<OpenMS::TargetedExperimentHelper::Peptide>& rhs) { return self.setPeptides(rhs); }, "rhs"_a)
        .def("getPeptides", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Peptide> & { return self.getPeptides(); }, nb::rv_policy::reference_internal)
        .def("hasPeptide", [](const OpenMS::TargetedExperiment& self, const OpenMS::String& ref) { return self.hasPeptide(ref); }, "ref"_a)
        .def("getPeptideByRef", [](const OpenMS::TargetedExperiment& self, const OpenMS::String& ref) -> const OpenMS::TargetedExperimentHelper::Peptide & { return self.getPeptideByRef(ref); }, "ref"_a, nb::rv_policy::reference_internal)
        .def("hasCompound", [](const OpenMS::TargetedExperiment& self, const OpenMS::String& ref) { return self.hasCompound(ref); }, "ref"_a)
        .def("getCompoundByRef", [](const OpenMS::TargetedExperiment& self, const OpenMS::String& ref) -> const OpenMS::TargetedExperimentHelper::Compound & { return self.getCompoundByRef(ref); }, "ref"_a, nb::rv_policy::reference_internal)
        .def("addPeptide", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperimentHelper::Peptide& rhs) { return self.addPeptide(rhs); }, "rhs"_a)
        .def("setTransitions", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::ReactionMonitoringTransition>& transitions) { return self.setTransitions(transitions); }, "transitions"_a)
        .def("setTransitions", [](OpenMS::TargetedExperiment& self, std::vector<OpenMS::ReactionMonitoringTransition>& transitions) { return self.setTransitions(transitions); }, "transitions"_a)
        .def("getTransitions", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::ReactionMonitoringTransition> & { return self.getTransitions(); }, nb::rv_policy::reference_internal)
        .def("addTransition", [](OpenMS::TargetedExperiment& self, const OpenMS::ReactionMonitoringTransition& transition) { return self.addTransition(transition); }, "transition"_a)
        .def("setIncludeTargets", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::IncludeExcludeTarget>& targets) { return self.setIncludeTargets(targets); }, "targets"_a)
        .def("getIncludeTargets", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::IncludeExcludeTarget> & { return self.getIncludeTargets(); }, nb::rv_policy::reference_internal)
        .def("addIncludeTarget", [](OpenMS::TargetedExperiment& self, const OpenMS::IncludeExcludeTarget& target) { return self.addIncludeTarget(target); }, "target"_a)
        .def("setExcludeTargets", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::IncludeExcludeTarget>& targets) { return self.setExcludeTargets(targets); }, "targets"_a)
        .def("getExcludeTargets", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::IncludeExcludeTarget> & { return self.getExcludeTargets(); }, nb::rv_policy::reference_internal)
        .def("addExcludeTarget", [](OpenMS::TargetedExperiment& self, const OpenMS::IncludeExcludeTarget& target) { return self.addExcludeTarget(target); }, "target"_a)
        .def("setSourceFiles", [](OpenMS::TargetedExperiment& self, const std::vector<OpenMS::SourceFile>& source_files) { return self.setSourceFiles(source_files); }, "source_files"_a)
        .def("getSourceFiles", [](const OpenMS::TargetedExperiment& self) -> const std::vector<OpenMS::SourceFile> & { return self.getSourceFiles(); }, nb::rv_policy::reference_internal)
        .def("addSourceFile", [](OpenMS::TargetedExperiment& self, const OpenMS::SourceFile& source_file) { return self.addSourceFile(source_file); }, "source_file"_a)
        .def("sortTransitionsByProductMZ", [](OpenMS::TargetedExperiment& self) { return self.sortTransitionsByProductMZ(); })
        .def("containsInvalidReferences", [](const OpenMS::TargetedExperiment& self) { return self.containsInvalidReferences(); })
        .def("__iadd__", [](OpenMS::TargetedExperiment& self, const OpenMS::TargetedExperiment& rhs) -> OpenMS::TargetedExperiment& { self += rhs; return self; }, "rhs"_a, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // TransformationModel_DataPoint
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationModel::DataPoint>(m, "TransformationModel_DataPoint",
        "Coordinate pair (with optional annotation) used for transformation models")
        .def(nb::init<>())
        .def(nb::init<double, double, const OpenMS::String&>(), "first"_a = 0.0, "second"_a = 0.0, "note"_a = "")
        .def(nb::init<const std::pair<double, double>&>(), "pair"_a)
        .def(nb::init<const OpenMS::TransformationModel::DataPoint&>())
        .def("__copy__", [](const OpenMS::TransformationModel::DataPoint& self) { return OpenMS::TransformationModel::DataPoint(self); })
        .def("__deepcopy__", [](const OpenMS::TransformationModel::DataPoint& self, nb::dict) { return OpenMS::TransformationModel::DataPoint(self); }, "memo"_a)
        .def_rw("first", &OpenMS::TransformationModel::DataPoint::first)
        .def_rw("second", &OpenMS::TransformationModel::DataPoint::second)
        .def_rw("note", &OpenMS::TransformationModel::DataPoint::note)
        .def(nb::self < nb::self)
        .def(nb::self == nb::self)
        ;
    // Alias for backwards compatibility with Cython API
    m.attr("TM_DataPoint") = m.attr("TransformationModel_DataPoint");

    // -----------------------------------------------------------------------
    // TransformationStatistics (TransformationDescription::TransformationStatistics)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationDescription::TransformationStatistics>(m, "TransformationStatistics",
        "Statistics for a transformation including percentile deviations")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TransformationDescription::TransformationStatistics &>())
        .def("__copy__", [](const OpenMS::TransformationDescription::TransformationStatistics& self) { return OpenMS::TransformationDescription::TransformationStatistics(self); })
        .def("__deepcopy__", [](const OpenMS::TransformationDescription::TransformationStatistics& self, nb::dict) { return OpenMS::TransformationDescription::TransformationStatistics(self); }, "memo"_a)
        .def_rw("percents", &OpenMS::TransformationDescription::TransformationStatistics::percents)
        .def_rw("xmin", &OpenMS::TransformationDescription::TransformationStatistics::xmin)
        .def_rw("xmax", &OpenMS::TransformationDescription::TransformationStatistics::xmax)
        .def_rw("ymin", &OpenMS::TransformationDescription::TransformationStatistics::ymin)
        .def_rw("ymax", &OpenMS::TransformationDescription::TransformationStatistics::ymax)
        .def_rw("percentiles_before", &OpenMS::TransformationDescription::TransformationStatistics::percentiles_before)
        .def_rw("percentiles_after", &OpenMS::TransformationDescription::TransformationStatistics::percentiles_after)
        ;

    // -----------------------------------------------------------------------
    // TransformationDescription
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationDescription>(m, "TransformationDescription", "Generic description of a coordinate transformation")
        .def(nb::init<>())
        .def(nb::init<std::vector<OpenMS::TransformationModel::DataPoint>>())
        .def(nb::init<const OpenMS::TransformationDescription &>())
        .def("__copy__", [](const OpenMS::TransformationDescription& self) { return OpenMS::TransformationDescription(self); })
        .def("__deepcopy__", [](const OpenMS::TransformationDescription& self, nb::dict) { return OpenMS::TransformationDescription(self); }, "memo"_a)
        .def("getModelType", [](const OpenMS::TransformationDescription& self) { return self.getModelType(); }, "Gets the type of the fitted model")
        .def("getModelParameters", [](const OpenMS::TransformationDescription& self) -> const OpenMS::Param & { return self.getModelParameters(); }, nb::rv_policy::reference_internal, "Returns the model parameters")
        .def("invert", [](OpenMS::TransformationDescription& self) { return self.invert(); }, "Computes an (approximate) inverse of the transformation")
        .def("getDeviations", [](const OpenMS::TransformationDescription& self, bool do_apply, bool do_sort) { std::vector<double> diffs; self.getDeviations(diffs, do_apply, do_sort); return diffs; }, "do_apply"_a = false, "do_sort"_a = true)
        .def("getStatistics", [](const OpenMS::TransformationDescription& self) { return self.getStatistics(); }, 
            R"doc(
Get the deviations between the data pairs
:param diffs: Output
:param do_apply: Get deviations after applying the model?
:param do_sort: Sort `diffs` before returning?
)doc")

        .def("getDataPoints", [](const OpenMS::TransformationDescription& self) {
            return self.getDataPoints();
        }, "Get the data points used for the transformation")

        .def("setDataPoints", [](OpenMS::TransformationDescription& self, const std::vector<std::pair<double, double>>& data) {
            OpenMS::TransformationDescription::DataPoints dp;
            for (const auto& p : data) {
                dp.push_back(std::make_pair(p.first, p.second));
            }
            self.setDataPoints(dp);
        }, "data"_a, "Set the data points for the transformation (from pairs)")

        .def("setDataPoints", [](OpenMS::TransformationDescription& self, const std::vector<OpenMS::TransformationModel::DataPoint>& data) {
            self.setDataPoints(data);
        }, "data"_a, "Set the data points for the transformation (from DataPoint objects)")

        .def("apply", [](const OpenMS::TransformationDescription& self, double value) {
            return self.apply(value);
        }, "value"_a, "Apply the transformation to a value")

        .def("fitModel", [](OpenMS::TransformationDescription& self, const OpenMS::String& model_type, const OpenMS::Param& params) {
            self.fitModel(model_type, params);
        }, "model_type"_a, "params"_a = OpenMS::Param(), "Fits a model to the data")

        .def_static("getModelTypes", [](nb::list result) {
            std::vector<OpenMS::String> types;
            OpenMS::TransformationDescription::getModelTypes(types);
            for (const auto& t : types) {
                result.append(nb::str(t.c_str()));
            }
        }, "result"_a, "Get available model types (fills list)")
        ;

    // -----------------------------------------------------------------------
    // TransformationModelLinear
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationModelLinear>(m, "TransformationModelLinear", 
        R"doc(
Linear model for transformations
TransformationModel
)doc")
        .def(nb::init<std::vector<OpenMS::TransformationModel::DataPoint>, OpenMS::Param>())
        .def("evaluate", [](const OpenMS::TransformationModelLinear& self, double value) { return self.evaluate(value); }, "value"_a)
        .def("invert", [](OpenMS::TransformationModelLinear& self) { return self.invert(); })
        .def("weightData", [](OpenMS::TransformationModelLinear& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.weightData(data);
            return data;
        }, "data"_a, "Weight the data by the given weight function, returns weighted data")
        .def("unWeightData", [](OpenMS::TransformationModelLinear& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.unWeightData(data);
            return data;
        }, "data"_a, "Unweight the data by the given weight function, returns unweighted data")
        .def("checkValidWeight", [](const OpenMS::TransformationModelLinear& self, const OpenMS::String& weight, const std::vector<OpenMS::String>& valid_weights) { return self.checkValidWeight(weight, valid_weights); }, "weight"_a, "valid_weights"_a, "Check for a valid weighting function string")
        .def("checkDatumRange", [](OpenMS::TransformationModelLinear& self, const double& datum, const double& datum_min, const double& datum_max) { return self.checkDatumRange(datum, datum_min, datum_max); }, "datum"_a, "datum_min"_a, "datum_max"_a, "Check that the datum is within the valid min and max bounds")
        .def("weightDatum", [](const OpenMS::TransformationModelLinear& self, const double& datum, const OpenMS::String& weight) { return self.weightDatum(datum, weight); }, "datum"_a, "weight"_a, "Weight the data according to the weighting function")
        .def("unWeightDatum", [](const OpenMS::TransformationModelLinear& self, const double& datum, const OpenMS::String& weight) { return self.unWeightDatum(datum, weight); }, "datum"_a, "weight"_a, "Apply the reverse of the weighting function to the data")
        .def("getValidXWeights", [](const OpenMS::TransformationModelLinear& self) { return self.getValidXWeights(); }, "Returns a list of valid x weight function stringss")
        .def("getValidYWeights", [](const OpenMS::TransformationModelLinear& self) { return self.getValidYWeights(); }, "Returns a list of valid y weight function strings")
        .def("getParameters", [](const OpenMS::TransformationModelLinear& self) -> const OpenMS::Param& { return self.getParameters(); }, nb::rv_policy::reference_internal, "Gets the (actual) parameters")

        .def_static("getDefaultParameters", [](OpenMS::Param& params) {
            OpenMS::TransformationModelLinear::getDefaultParameters(params);
        }, "params"_a, "Get default parameters")
        ;

    // -----------------------------------------------------------------------
    // XQuestScores
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::XQuestScores>(m, "XQuestScores", 
        R"doc(
An implementation of the scores for cross-link identification from
the xQuest algorithm (O. Rinner et al., 2008, "Identification of
cross-linked peptides from large sequence databases")
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::XQuestScores &>())
        .def("__copy__", [](const OpenMS::XQuestScores& self) { return OpenMS::XQuestScores(self); })
        .def("__deepcopy__", [](const OpenMS::XQuestScores& self, nb::dict) { return OpenMS::XQuestScores(self); }, "memo"_a)
        .def_static("preScore", [](size_t matched_alpha, size_t ions_alpha, size_t matched_beta, size_t ions_beta) { return OpenMS::XQuestScores::preScore(matched_alpha, ions_alpha, matched_beta, ions_beta); }, "matched_alpha"_a, "ions_alpha"_a, "matched_beta"_a, "ions_beta"_a)
        .def_static("preScore", [](size_t matched_alpha, size_t ions_alpha) { return OpenMS::XQuestScores::preScore(matched_alpha, ions_alpha); }, "matched_alpha"_a, "ions_alpha"_a, 
            R"doc(
Compute a simple and fast to compute pre-score for a cross-link spectrum match
:param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
:param ions_alpha: Number of theoretical ions from the alpha peptide
:param matched_beta: Number of experimental peaks matched to theoretical linear ions from the beta peptide
:param ions_beta: Number of theoretical ions from the beta peptide
)doc")
        .def_static("matchOddsScore", [](const OpenMS::MSSpectrum& theoretical_spec, size_t matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum, size_t n_charges) { return OpenMS::XQuestScores::matchOddsScore(theoretical_spec, matched_size, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, is_xlink_spectrum, n_charges); }, "theoretical_spec"_a, "matched_size"_a, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, "is_xlink_spectrum"_a, "n_charges"_a, 
            R"doc(
Compute a simple and fast to compute pre-score for a mono-link spectrum match
:param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
:param ions_alpha: Number of theoretical ions from the alpha peptide
)doc")
        .def_static("logOccupancyProb", [](const OpenMS::MSSpectrum& theoretical_spec, size_t matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm) { return OpenMS::XQuestScores::logOccupancyProb(theoretical_spec, matched_size, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm); }, "theoretical_spec"_a, "matched_size"_a, "fragment_mass_tolerance"_a, "fragment_mass_tolerance_unit_ppm"_a, 
            R"doc(
Compute the match-odds score, a score based on the probability of getting the given number of matched peaks by chance
:param theoretical_spec: Theoretical spectrum, sorted by position
:param matched_size: Alignment between the theoretical and the experimental spectra
:param fragment_mass_tolerance: Fragment mass tolerance of the alignment
:param fragment_mass_tolerance_unit_ppm: Fragment mass tolerance unit of the alignment, true = ppm, false = Da
:param is_xlink_spectrum: Type of cross-link, true = cross-link, false = mono-link
:param n_charges: Number of considered charges in the theoretical spectrum
)doc")
        .def_static("weightedTICScoreXQuest", [](size_t alpha_size, size_t beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link) { return OpenMS::XQuestScores::weightedTICScoreXQuest(alpha_size, beta_size, intsum_alpha, intsum_beta, total_current, type_is_cross_link); }, "alpha_size"_a, "beta_size"_a, "intsum_alpha"_a, "intsum_beta"_a, "total_current"_a, "type_is_cross_link"_a, 
            R"doc(
Compute the logOccupancyProb score, similar to the match_odds, a score based on the probability of getting the given number of matched peaks by chance
:param theoretical_spec: Theoretical spectrum, sorted by position
:param matched_size: Number of matched peaks between experimental and theoretical spectra
:param fragment_mass_tolerance: The tolerance of the alignment
:param fragment_mass_tolerance_unit: The tolerance unit of the alignment, true = ppm, false = Da
)doc")
        .def_static("weightedTICScore", [](size_t alpha_size, size_t beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link) { return OpenMS::XQuestScores::weightedTICScore(alpha_size, beta_size, intsum_alpha, intsum_beta, total_current, type_is_cross_link); }, "alpha_size"_a, "beta_size"_a, "intsum_alpha"_a, "intsum_beta"_a, "total_current"_a, "type_is_cross_link"_a)
        .def_static("xCorrelation", [](const OpenMS::MSSpectrum& spec1, const OpenMS::MSSpectrum& spec2, int maxshift, double tolerance) { return OpenMS::XQuestScores::xCorrelation(spec1, spec2, maxshift, tolerance); }, "spec1"_a, "spec2"_a, "maxshift"_a, "tolerance"_a)
        .def_static("xCorrelationPrescore", [](const OpenMS::MSSpectrum& spec1, const OpenMS::MSSpectrum& spec2, double tolerance) { return OpenMS::XQuestScores::xCorrelationPrescore(spec1, spec2, tolerance); }, "spec1"_a, "spec2"_a, "tolerance"_a)
        .def_static("matchedCurrentChain", [](const std::vector<std::pair<size_t, size_t>>& matched_spec_linear, const std::vector<std::pair<size_t, size_t>>& matched_spec_xlinks, const OpenMS::MSSpectrum& spectrum_linear_peaks, const OpenMS::MSSpectrum& spectrum_xlink_peaks) { return OpenMS::XQuestScores::matchedCurrentChain(matched_spec_linear, matched_spec_xlinks, spectrum_linear_peaks, spectrum_xlink_peaks); }, "matched_spec_linear"_a, "matched_spec_xlinks"_a, "spectrum_linear_peaks"_a, "spectrum_xlink_peaks"_a, "Computes sum of peak intensities of matched peaks for either alpha or beta peptide")
        .def_static("totalMatchedCurrent", [](const std::vector<std::pair<size_t, size_t>>& matched_spec_linear_alpha, const std::vector<std::pair<size_t, size_t>>& matched_spec_linear_beta, const std::vector<std::pair<size_t, size_t>>& matched_spec_xlinks_alpha, const std::vector<std::pair<size_t, size_t>>& matched_spec_xlinks_beta, const OpenMS::MSSpectrum& spectrum_linear_peaks, const OpenMS::MSSpectrum& spectrum_xlink_peaks) { return OpenMS::XQuestScores::totalMatchedCurrent(matched_spec_linear_alpha, matched_spec_linear_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta, spectrum_linear_peaks, spectrum_xlink_peaks); }, "matched_spec_linear_alpha"_a, "matched_spec_linear_beta"_a, "matched_spec_xlinks_alpha"_a, "matched_spec_xlinks_beta"_a, "spectrum_linear_peaks"_a, "spectrum_xlink_peaks"_a, "Computes sum of peak intensities of all matched peaks")
        ;


    using OSBDA = OpenSwath::OSBinaryDataArray;
    // -----------------------------------------------------------------------
    // OSBinaryDataArray
    // -----------------------------------------------------------------------
    nb::class_<OSBDA>(m, "OSBinaryDataArray", "OpenMS class OSBinaryDataArray")
        .def(nb::init<>())
        .def(nb::init<const OSBDA&>())
        .def("__copy__", [](const OSBDA& self) { return OSBDA(self); })
        .def("__deepcopy__", [](const OSBDA& self, nb::dict) { return OSBDA(self); }, "memo"_a)
        .def_rw("data", &OSBDA::data)
        .def_rw("description", &OSBDA::description)
        .def("get_data", [](const OSBDA& self) {
            return self.data;
        }, "Access to a copy of the underlying data")
        .def("get_data_view", [](nb::object self_obj) -> nb::object {
            auto& self = nb::cast<OSBDA&>(self_obj);
            double* data_ptr = self.data.empty() ? nullptr : self.data.data();
            size_t shape[] = {self.data.size()};
            return nb::ndarray<nb::numpy, double>(data_ptr, 1, shape, self_obj).cast();
        }, "Access to the underlying data using a writable view (empty array if empty)")
        ;


    // -----------------------------------------------------------------------
    // OSSpectrumMeta
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::OSSpectrumMeta>(m, "OSSpectrumMeta", "Identifying information for a spectrum")
        .def(nb::init<>())
        .def("__copy__", [](const OpenSwath::OSSpectrumMeta& self) { return OpenSwath::OSSpectrumMeta(self); })
        .def("__deepcopy__", [](const OpenSwath::OSSpectrumMeta& self, nb::dict) { return OpenSwath::OSSpectrumMeta(self); }, "memo"_a)
        .def_rw("index", &OpenSwath::OSSpectrumMeta::index)
        .def_rw("id", &OpenSwath::OSSpectrumMeta::id)
        .def_rw("RT", &OpenSwath::OSSpectrumMeta::RT)
        .def_rw("ms_level", &OpenSwath::OSSpectrumMeta::ms_level)
        ;

    // -----------------------------------------------------------------------
    // ChromExtractParams (OSW_ChromExtractParams)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromExtractParams>(m, "OSW_ChromExtractParams",
        "Parameters for chromatogram extraction in OpenSWATH workflow")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::ChromExtractParams& self) { return OpenMS::ChromExtractParams(self); })
        .def("__deepcopy__", [](const OpenMS::ChromExtractParams& self, nb::dict) { return OpenMS::ChromExtractParams(self); }, "memo"_a)
        .def_rw("min_upper_edge_dist", &OpenMS::ChromExtractParams::min_upper_edge_dist)
        .def_rw("mz_extraction_window", &OpenMS::ChromExtractParams::mz_extraction_window)
        .def_rw("im_extraction_window", &OpenMS::ChromExtractParams::im_extraction_window)
        .def_rw("ppm", &OpenMS::ChromExtractParams::ppm)
        .def_rw("extraction_function", &OpenMS::ChromExtractParams::extraction_function)
        .def_rw("rt_extraction_window", &OpenMS::ChromExtractParams::rt_extraction_window)
        .def_rw("extra_rt_extract", &OpenMS::ChromExtractParams::extra_rt_extract)
        ;

    // -----------------------------------------------------------------------
    // OSChromatogramMeta
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::OSChromatogramMeta>(m, "OSChromatogramMeta", "Identifying information for a chromatogram")
        .def(nb::init<>())
        .def("__copy__", [](const OpenSwath::OSChromatogramMeta& self) { return OpenSwath::OSChromatogramMeta(self); })
        .def("__deepcopy__", [](const OpenSwath::OSChromatogramMeta& self, nb::dict) { return OpenSwath::OSChromatogramMeta(self); }, "memo"_a)
        .def_rw("index", &OpenSwath::OSChromatogramMeta::index)
        .def_rw("id", &OpenSwath::OSChromatogramMeta::id)
        ;

    using OSSpec = OpenSwath::OSSpectrum;
    using BDAPtr = std::shared_ptr<OpenSwath::BinaryDataArray>;
    // -----------------------------------------------------------------------
    // OSSpectrum
    // -----------------------------------------------------------------------
    nb::class_<OSSpec>(m, "OSSpectrum", "OpenMS class OSSpectrum")
        .def(nb::init<>())
        .def(nb::init<const OSSpec&>())
        .def("__copy__", [](const OSSpec& self) { return OSSpec(self); })
        .def("__deepcopy__", [](const OSSpec& self, nb::dict) { return OSSpec(self); }, "memo"_a)
        .def("getMZArray", &OSSpec::getMZArray)
        .def("getIntensityArray", &OSSpec::getIntensityArray)
        .def("getDriftTimeArray", &OSSpec::getDriftTimeArray)
        .def("setMZArray", &OSSpec::setMZArray, "data"_a)
        .def("setIntensityArray", &OSSpec::setIntensityArray, "data"_a)
        .def("setDriftTimeArray", &OSSpec::setDriftTimeArray, "data"_a)
        .def("set_mz_array", [](OSSpec& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA>();
            arr->data = std::move(data);
            self.setMZArray(arr);
        }, "data"_a, "Set m/z array from list")
        .def("set_intensity_array", [](OSSpec& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")
        .def("get_mz_array", [](const OSSpec& self) {
            auto arr = self.getMZArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get m/z array as list")
        .def("get_mz_array_view", [](OSSpec& self) -> nb::object {
            auto mz_arr = self.getMZArray();
            if (!mz_arr || mz_arr->data.empty()) {
                size_t shape[] = {0};
                return nb::ndarray<nb::numpy, double>(nullptr, 1, shape, nb::none()).cast();
            }
            nb::object owner = nb::cast(mz_arr);
            auto& data = mz_arr->data;
            size_t shape[] = {data.size()};
            return nb::ndarray<nb::numpy, double>(data.data(), 1, shape, owner).cast();
        }, "Get m/z array as writable view (empty array if no data)")
        .def("get_intensity_array", [](const OSSpec& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array as list")
        .def("get_intensity_array_view", [](OSSpec& self) -> nb::object {
            auto int_arr = self.getIntensityArray();
            if (!int_arr || int_arr->data.empty()) {
                size_t shape[] = {0};
                return nb::ndarray<nb::numpy, double>(nullptr, 1, shape, nb::none()).cast();
            }
            nb::object owner = nb::cast(int_arr);
            auto& data = int_arr->data;
            size_t shape[] = {data.size()};
            return nb::ndarray<nb::numpy, double>(data.data(), 1, shape, owner).cast();
        }, "Get intensity array as writable view (empty array if no data)")
        .def("get_drift_time_array", [](const OSSpec& self) -> nb::object {
            auto arr = self.getDriftTimeArray();
            if (!arr) return nb::none();
            return nb::cast(arr->data);
        }, "Get drift time array or None")
        .def("get_drift_time_array_view", [](OSSpec& self) -> nb::object {
            auto arr = self.getDriftTimeArray();
            if (!arr || arr->data.empty()) {
                size_t shape[] = {0};
                return nb::ndarray<nb::numpy, double>(nullptr, 1, shape, nb::none()).cast();
            }
            nb::object owner = nb::cast(arr);
            auto& data = arr->data;
            size_t shape[] = {data.size()};
            return nb::ndarray<nb::numpy, double>(data.data(), 1, shape, owner).cast();
        }, "Get drift time array as writable view (empty array if no data)")
        .def("get_data_arrays", [](OSSpec& self) {
            auto& arrays = self.getDataArrays();
            std::vector<std::shared_ptr<OSBDA>> result;
            for (auto& a : arrays) result.push_back(a);
            return result;
        }, "Get all data arrays")
        .def("set_data_arrays", [](OSSpec& self, std::vector<std::shared_ptr<OSBDA>> arrays) {
            std::vector<BDAPtr> ptrs;
            for (auto& a : arrays) ptrs.push_back(a);
            self.setDataArrays(ptrs);
        }, "arrays"_a, "Set all data arrays")
        ;


    using OSChrom = OpenSwath::OSChromatogram;
    using OSBDA_C = OpenSwath::OSBinaryDataArray;
    using BDAPtr_C = std::shared_ptr<OpenSwath::BinaryDataArray>;
    // -----------------------------------------------------------------------
    // OSChromatogram
    // -----------------------------------------------------------------------
    nb::class_<OSChrom>(m, "OSChromatogram", "OpenMS class OSChromatogram")
        .def(nb::init<>())
        .def(nb::init<const OSChrom&>())
        .def("__copy__", [](const OSChrom& self) { return OSChrom(self); })
        .def("__deepcopy__", [](const OSChrom& self, nb::dict) { return OSChrom(self); }, "memo"_a)
        .def("getTimeArray", &OSChrom::getTimeArray)
        .def("getIntensityArray", &OSChrom::getIntensityArray)
        .def("setTimeArray", &OSChrom::setTimeArray, "data"_a)
        .def("setIntensityArray", &OSChrom::setIntensityArray, "data"_a)
        .def("set_time_array", [](OSChrom& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA_C>();
            arr->data = std::move(data);
            self.setTimeArray(arr);
        }, "data"_a, "Set time array from list")
        .def("set_intensity_array", [](OSChrom& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA_C>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")
        .def("get_time_array", [](OSChrom& self) {
            auto arr = self.getTimeArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get time array as list")
        .def("get_intensity_array", [](OSChrom& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array as list")
        .def("get_time_array_view", [](OSChrom& self) -> nb::object {
            auto time_arr = self.getTimeArray();
            if (!time_arr || time_arr->data.empty()) {
                size_t shape[] = {0};
                return nb::ndarray<nb::numpy, double>(nullptr, 1, shape, nb::none()).cast();
            }
            nb::object owner = nb::cast(time_arr);
            auto& data = time_arr->data;
            size_t shape[] = {data.size()};
            return nb::ndarray<nb::numpy, double>(data.data(), 1, shape, owner).cast();
        }, "Get time array as writable view (empty array if no data)")
        .def("get_intensity_array_view", [](OSChrom& self) -> nb::object {
            auto int_arr = self.getIntensityArray();
            if (!int_arr || int_arr->data.empty()) {
                size_t shape[] = {0};
                return nb::ndarray<nb::numpy, double>(nullptr, 1, shape, nb::none()).cast();
            }
            nb::object owner = nb::cast(int_arr);
            auto& data = int_arr->data;
            size_t shape[] = {data.size()};
            return nb::ndarray<nb::numpy, double>(data.data(), 1, shape, owner).cast();
        }, "Get intensity array as writable view (empty array if no data)")
        .def("get_data_arrays", [](OSChrom& self) {
            auto& arrays = self.getDataArrays();
            std::vector<std::shared_ptr<OSBDA_C>> result;
            for (auto& a : arrays) result.push_back(a);
            return result;
        }, "Get all data arrays")
        .def("set_data_arrays", [](OSChrom& self, std::vector<std::shared_ptr<OSBDA_C>> arrays) {
            std::vector<BDAPtr_C> ptrs;
            for (auto& a : arrays) ptrs.push_back(a);
            self.setDataArrays(ptrs);
        }, "arrays"_a, "Set all data arrays")
        ;


    using ICI = OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation;
    // -----------------------------------------------------------------------
    // IsobaricChannelInformation
    // -----------------------------------------------------------------------
    nb::class_<ICI>(m, "IsobaricChannelInformation", "OpenMS class IsobaricChannelInformation")
        .def("__init__", [](ICI* self, const std::string& name, int id, const std::string& description, double center, std::vector<int> affected_channels) {
            new (self) ICI(name, id, description, center, affected_channels);
        }, "name"_a, "id"_a, "description"_a, "center"_a, "affected_channels"_a)
        .def(nb::init<const ICI&>())
        .def("__copy__", [](const ICI& self) { return ICI(self); })
        .def("__deepcopy__", [](const ICI& self, nb::dict) { return ICI(self); }, "memo"_a)
        .def_rw("name", &ICI::name)
        .def_rw("id", &ICI::id)
        .def_rw("description", &ICI::description)
        .def_rw("center", &ICI::center)
        .def_rw("affected_channels", &ICI::affected_channels)
        ;


    // -----------------------------------------------------------------------
    // TransformationModelBSpline
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationModelBSpline>(m, "TransformationModelBSpline", "B-spline model for transformations")
        .def(nb::init<const std::vector<OpenMS::TransformationModel::DataPoint>&, const OpenMS::Param&>(), "data"_a, "params"_a)
        .def("evaluate", [](const OpenMS::TransformationModelBSpline& self, double value) { return self.evaluate(value); }, "value"_a)
        .def("getParameters", [](const OpenMS::TransformationModelBSpline& self) -> const OpenMS::Param& { return self.getParameters(); }, nb::rv_policy::reference_internal, "Gets the (actual) parameters")
        .def("weightData", [](OpenMS::TransformationModelBSpline& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.weightData(data);
            return data;
        }, "data"_a, "Weight the data by the given weight function, returns weighted data")
        .def("unWeightData", [](OpenMS::TransformationModelBSpline& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.unWeightData(data);
            return data;
        }, "data"_a, "Unweight the data by the given weight function, returns unweighted data")
        .def("checkValidWeight", [](const OpenMS::TransformationModelBSpline& self, const OpenMS::String& weight, const std::vector<OpenMS::String>& valid_weights) { return self.checkValidWeight(weight, valid_weights); }, "weight"_a, "valid_weights"_a, "Check for a valid weighting function string")
        .def("checkDatumRange", [](OpenMS::TransformationModelBSpline& self, const double& datum, const double& datum_min, const double& datum_max) { return self.checkDatumRange(datum, datum_min, datum_max); }, "datum"_a, "datum_min"_a, "datum_max"_a, "Check that the datum is within the valid min and max bounds")
        .def("weightDatum", [](const OpenMS::TransformationModelBSpline& self, const double& datum, const OpenMS::String& weight) { return self.weightDatum(datum, weight); }, "datum"_a, "weight"_a, "Weight the data according to the weighting function")
        .def("unWeightDatum", [](const OpenMS::TransformationModelBSpline& self, const double& datum, const OpenMS::String& weight) { return self.unWeightDatum(datum, weight); }, "datum"_a, "weight"_a, "Apply the reverse of the weighting function to the data")
        .def("getValidXWeights", [](const OpenMS::TransformationModelBSpline& self) { return self.getValidXWeights(); }, "Returns a list of valid x weight function strings")
        .def("getValidYWeights", [](const OpenMS::TransformationModelBSpline& self) { return self.getValidYWeights(); }, "Returns a list of valid y weight function strings")
        .def_static("getDefaultParameters", [](OpenMS::Param& params) {
            OpenMS::TransformationModelBSpline::getDefaultParameters(params);
        }, "params"_a, "Get default parameters")
        ;

    // -----------------------------------------------------------------------
    // TransformationModelLowess
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationModelLowess>(m, "TransformationModelLowess", "Lowess model for transformations")
        .def(nb::init<const std::vector<OpenMS::TransformationModel::DataPoint>&, const OpenMS::Param&>(), "data"_a, "params"_a)
        .def("evaluate", [](const OpenMS::TransformationModelLowess& self, double value) { return self.evaluate(value); }, "value"_a)
        .def("getParameters", [](const OpenMS::TransformationModelLowess& self) -> const OpenMS::Param& { return self.getParameters(); }, nb::rv_policy::reference_internal, "Gets the (actual) parameters")
        .def("weightData", [](OpenMS::TransformationModelLowess& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.weightData(data);
            return data;
        }, "data"_a, "Weight the data by the given weight function, returns weighted data")
        .def("unWeightData", [](OpenMS::TransformationModelLowess& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.unWeightData(data);
            return data;
        }, "data"_a, "Unweight the data by the given weight function, returns unweighted data")
        .def("checkValidWeight", [](const OpenMS::TransformationModelLowess& self, const OpenMS::String& weight, const std::vector<OpenMS::String>& valid_weights) { return self.checkValidWeight(weight, valid_weights); }, "weight"_a, "valid_weights"_a, "Check for a valid weighting function string")
        .def("checkDatumRange", [](OpenMS::TransformationModelLowess& self, const double& datum, const double& datum_min, const double& datum_max) { return self.checkDatumRange(datum, datum_min, datum_max); }, "datum"_a, "datum_min"_a, "datum_max"_a, "Check that the datum is within the valid min and max bounds")
        .def("weightDatum", [](const OpenMS::TransformationModelLowess& self, const double& datum, const OpenMS::String& weight) { return self.weightDatum(datum, weight); }, "datum"_a, "weight"_a, "Weight the data according to the weighting function")
        .def("unWeightDatum", [](const OpenMS::TransformationModelLowess& self, const double& datum, const OpenMS::String& weight) { return self.unWeightDatum(datum, weight); }, "datum"_a, "weight"_a, "Apply the reverse of the weighting function to the data")
        .def("getValidXWeights", [](const OpenMS::TransformationModelLowess& self) { return self.getValidXWeights(); }, "Returns a list of valid x weight function strings")
        .def("getValidYWeights", [](const OpenMS::TransformationModelLowess& self) { return self.getValidYWeights(); }, "Returns a list of valid y weight function strings")
        .def_static("getDefaultParameters", [](OpenMS::Param& params) {
            OpenMS::TransformationModelLowess::getDefaultParameters(params);
        }, "params"_a, "Get default parameters")
        ;


    // -----------------------------------------------------------------------
    // IsobaricNormalizer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsobaricNormalizer>(m, "IsobaricNormalizer", "OpenMS class IsobaricNormalizer")
        .def("__copy__", [](const OpenMS::IsobaricNormalizer& self) { return OpenMS::IsobaricNormalizer(self); })
        .def("__deepcopy__", [](const OpenMS::IsobaricNormalizer& self, nb::dict) { return OpenMS::IsobaricNormalizer(self); }, "memo"_a)
        .def(nb::init<const OpenMS::ItraqFourPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::ItraqEightPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTSixPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTTenPlexQuantitationMethod*>(), "quant_method"_a)
        .def("normalize", &OpenMS::IsobaricNormalizer::normalize, "consensus_map"_a, "Normalize consensus map")
        ;


    // -----------------------------------------------------------------------
    // SpectrumAccessOpenMS
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumAccessOpenMS>(m, "SpectrumAccessOpenMS", "OpenMS class SpectrumAccessOpenMS")
        .def("__copy__", [](const OpenMS::SpectrumAccessOpenMS& self) { return OpenMS::SpectrumAccessOpenMS(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumAccessOpenMS& self, nb::dict) { return OpenMS::SpectrumAccessOpenMS(self); }, "memo"_a)
        .def(nb::init<std::shared_ptr<OpenMS::MSExperiment>>(), "ms_experiment"_a)
        .def("getNrSpectra", [](const OpenMS::SpectrumAccessOpenMS& self) { return self.getNrSpectra(); }, "Get number of spectra")
        .def("getNrChromatograms", [](const OpenMS::SpectrumAccessOpenMS& self) { return self.getNrChromatograms(); }, "Get number of chromatograms")
        .def("getSpectrumById", [](OpenMS::SpectrumAccessOpenMS& self, int id) { return self.getSpectrumById(id); }, "id"_a, "Get spectrum by index")
        .def("getChromatogramById", [](OpenMS::SpectrumAccessOpenMS& self, int id) { return self.getChromatogramById(id); }, "id"_a, "Get chromatogram by index")
        .def("getChromatogramNativeID", [](const OpenMS::SpectrumAccessOpenMS& self, int id) { return self.getChromatogramNativeID(id); }, "id"_a, "Returns the native ID of the chromatogram")
        .def("getSpectraByRT", [](const OpenMS::SpectrumAccessOpenMS& self, double RT, double deltaRT) { return self.getSpectraByRT(RT, deltaRT); }, "RT"_a, "deltaRT"_a, "Returns spectra indices within RT range")
        ;


    // -----------------------------------------------------------------------
    // SpectrumAccessOpenMSInMemory
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumAccessOpenMSInMemory>(m, "SpectrumAccessOpenMSInMemory", "OpenMS class SpectrumAccessOpenMSInMemory")
        .def("__copy__", [](const OpenMS::SpectrumAccessOpenMSInMemory& self) { return OpenMS::SpectrumAccessOpenMSInMemory(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumAccessOpenMSInMemory& self, nb::dict) { return OpenMS::SpectrumAccessOpenMSInMemory(self); }, "memo"_a)
        .def("__init__", [](OpenMS::SpectrumAccessOpenMSInMemory* self, OpenMS::SpectrumAccessOpenMS& other) {
            new (self) OpenMS::SpectrumAccessOpenMSInMemory(other);
        }, "other"_a)
        .def("getNrSpectra", [](const OpenMS::SpectrumAccessOpenMSInMemory& self) { return self.getNrSpectra(); })
        .def("getNrChromatograms", [](const OpenMS::SpectrumAccessOpenMSInMemory& self) { return self.getNrChromatograms(); })
        .def("getSpectrumById", [](OpenMS::SpectrumAccessOpenMSInMemory& self, int id) { return self.getSpectrumById(id); }, "id"_a, "Get spectrum by index")
        .def("getChromatogramById", [](OpenMS::SpectrumAccessOpenMSInMemory& self, int id) { return self.getChromatogramById(id); }, "id"_a, "Get chromatogram by index")
        .def("getChromatogramNativeID", [](const OpenMS::SpectrumAccessOpenMSInMemory& self, int id) { return self.getChromatogramNativeID(id); }, "id"_a, "Returns the native ID of the chromatogram")
        .def("getSpectraByRT", [](const OpenMS::SpectrumAccessOpenMSInMemory& self, double RT, double deltaRT) { return self.getSpectraByRT(RT, deltaRT); }, "RT"_a, "deltaRT"_a, "Returns spectra indices within RT range")
        ;


    // -----------------------------------------------------------------------
    // SwathMap
    // -----------------------------------------------------------------------
    nb::class_<OpenSwath::SwathMap>(m, "SwathMap", "OpenMS class SwathMap")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::SwathMap&>())
        .def("__copy__", [](const OpenSwath::SwathMap& self) { return OpenSwath::SwathMap(self); })
        .def("__deepcopy__", [](const OpenSwath::SwathMap& self, nb::dict) { return OpenSwath::SwathMap(self); }, "memo"_a)
        .def(nb::init<double, double, double, bool>(), "mz_start"_a, "mz_end"_a, "mz_center"_a, "is_ms1"_a)
        .def_rw("lower", &OpenSwath::SwathMap::lower)
        .def_rw("upper", &OpenSwath::SwathMap::upper)
        .def_rw("center", &OpenSwath::SwathMap::center)
        .def_rw("ms1", &OpenSwath::SwathMap::ms1)
        .def("setSpectrumPtr", [](OpenSwath::SwathMap& self, std::shared_ptr<OpenMS::SpectrumAccessOpenMS> sa) {
            self.sptr = sa;
        }, "spectrum_access"_a, "Set spectrum access pointer (SpectrumAccessOpenMS)")
        .def("setSpectrumPtr", [](OpenSwath::SwathMap& self, std::shared_ptr<OpenMS::SpectrumAccessOpenMSInMemory> sa) {
            self.sptr = sa;
        }, "spectrum_access"_a, "Set spectrum access pointer (SpectrumAccessOpenMSInMemory)")
        .def("getSpectrumPtr", [](OpenSwath::SwathMap& self) -> nb::object {
            if (!self.sptr) return nb::none();
            // Try dynamic_pointer_cast to known types
            auto sa = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOpenMS>(self.sptr);
            if (sa) return nb::cast(sa);
            auto inmem = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOpenMSInMemory>(self.sptr);
            if (inmem) return nb::cast(inmem);
            return nb::none();
        }, "Get spectrum access pointer")
        ;


    // OpenSwathScoring
    // -----------------------------------------------------------------------
    // OpenSwathScoring
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OpenSwathScoring>(m, "OpenSwathScoring", "OpenMS class OpenSwathScoring")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::OpenSwathScoring& self) { return OpenMS::OpenSwathScoring(self); })
        .def("__deepcopy__", [](const OpenMS::OpenSwathScoring& self, nb::dict) { return OpenMS::OpenSwathScoring(self); }, "memo"_a)
        .def("initialize", &OpenMS::OpenSwathScoring::initialize,
            "rt_normalization_factor"_a, "add_up_spectra"_a,
            "spacing_for_spectra_resampling"_a, "merge_spectra_by_peak_width_fraction"_a,
            "drift_extra"_a, "su"_a, "spectrum_addition_method"_a,
            "spectrum_merge_method_type"_a, "use_ms1_ion_mobility"_a,
            "apply_im_peak_picking"_a)
        ;


    // -----------------------------------------------------------------------
    // TransformationModelInterpolated
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TransformationModelInterpolated>(m, "TransformationModelInterpolated",
        R"doc(
Interpolation model for transformations

Between the data points, the interpolation uses the neighboring points.
Supports linear, cspline, and akima interpolation methods.
TransformationModel
)doc")
        .def(nb::init<const std::vector<OpenMS::TransformationModel::DataPoint>&, const OpenMS::Param&>(), "data"_a, "params"_a)
        .def("evaluate", [](const OpenMS::TransformationModelInterpolated& self, double value) { return self.evaluate(value); }, "value"_a)
        .def("getParameters", [](const OpenMS::TransformationModelInterpolated& self) -> const OpenMS::Param& { return self.getParameters(); }, nb::rv_policy::reference_internal, "Gets the (actual) parameters")
        .def("weightData", [](OpenMS::TransformationModelInterpolated& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.weightData(data);
            return data;
        }, "data"_a, "Weight the data by the given weight function, returns weighted data")
        .def("unWeightData", [](OpenMS::TransformationModelInterpolated& self, std::vector<OpenMS::TransformationModel::DataPoint> data) {
            self.unWeightData(data);
            return data;
        }, "data"_a, "Unweight the data by the given weight function, returns unweighted data")
        .def("checkValidWeight", [](const OpenMS::TransformationModelInterpolated& self, const OpenMS::String& weight, const std::vector<OpenMS::String>& valid_weights) { return self.checkValidWeight(weight, valid_weights); }, "weight"_a, "valid_weights"_a, "Check for a valid weighting function string")
        .def("checkDatumRange", [](OpenMS::TransformationModelInterpolated& self, const double& datum, const double& datum_min, const double& datum_max) { return self.checkDatumRange(datum, datum_min, datum_max); }, "datum"_a, "datum_min"_a, "datum_max"_a, "Check that the datum is within the valid min and max bounds")
        .def("weightDatum", [](const OpenMS::TransformationModelInterpolated& self, const double& datum, const OpenMS::String& weight) { return self.weightDatum(datum, weight); }, "datum"_a, "weight"_a, "Weight the data according to the weighting function")
        .def("unWeightDatum", [](const OpenMS::TransformationModelInterpolated& self, const double& datum, const OpenMS::String& weight) { return self.unWeightDatum(datum, weight); }, "datum"_a, "weight"_a, "Apply the reverse of the weighting function to the data")
        .def("getValidXWeights", [](const OpenMS::TransformationModelInterpolated& self) { return self.getValidXWeights(); }, "Returns a list of valid x weight function strings")
        .def("getValidYWeights", [](const OpenMS::TransformationModelInterpolated& self) { return self.getValidYWeights(); }, "Returns a list of valid y weight function strings")
        .def_static("getDefaultParameters", []() {
            OpenMS::Param params;
            OpenMS::TransformationModelInterpolated::getDefaultParameters(params);
            return params;
        }, "Get default parameters")
        ;

    // SpectrumAccessOpenMSCached, SpectrumAccessQuadMZTransforming: cannot bind
    // in analysis because they inherit from ISpectrumAccess/CachedmzML which
    // are not bound as nanobind base classes. Mark as xfail in tests.

    // -----------------------------------------------------------------------
    // ILPDCWrapper
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ILPDCWrapper>(m, "ILPDCWrapper",
        "Integer linear programming wrapper for decharging")
        .def(nb::init<>())
        .def("compute", [](const OpenMS::ILPDCWrapper& self, const OpenMS::FeatureMap& fm, std::vector<OpenMS::ChargePair>& pairs, OpenMS::Size verbose_level) {
            return self.compute(fm, pairs, verbose_level);
        }, "fm"_a, "pairs"_a, "verbose_level"_a, "Compute optimal charge pairs using ILP")
        ;

}