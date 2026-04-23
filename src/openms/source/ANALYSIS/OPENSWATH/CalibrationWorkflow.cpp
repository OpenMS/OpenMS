// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>
#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  CalibrationWorkflow::CalibrationWorkflow() :
    DefaultParamHandler("CalibrationWorkflow"),
    ProgressLogger()
  {
    // iRT peptide sampling parameters
    defaults_.setValue("auto_irt:enabled", "true", "Whether to sample iRTs on‐the‐fly (true) from the input targeted transition file (instead of passing specific iRT files). This may be useful if standard iRTs (Biognosys iRT kit) were not spiked-in. If set to false, and no additional iRT files are provided, and no transformation is provided, then no calibration is performed.");
    defaults_.setValidStrings("auto_irt:enabled", {"true", "false"});
    
    defaults_.setValue("auto_irt:irt_bins", 100, "Number of RT bins for linear iRT sampling");
    defaults_.setMinInt("auto_irt:irt_bins", 1);
    defaults_.setMaxInt("auto_irt:irt_bins", 10000);
    
    defaults_.setValue("auto_irt:irt_peptides_per_bin", 5, "Peptides sampled per bin for linear iRT");
    defaults_.setMinInt("auto_irt:irt_peptides_per_bin", 1);
    // defaults_.setMaxInt("auto_irt:irt_peptides_per_bin", 1000);
    
    defaults_.setValue("auto_irt:irt_seed", 5489, "RNG seed for reproducible sampling (0 = non-deterministic)");
    defaults_.setMinInt("auto_irt:irt_seed", 0);
    
    defaults_.setValue("auto_irt:irt_bins_nonlinear", 2000, "Number of RT bins for nonlinear iRT sampling");
    defaults_.setMinInt("auto_irt:irt_bins_nonlinear", 1);
    defaults_.setMaxInt("auto_irt:irt_bins_nonlinear", 10000);
    
    defaults_.setValue("auto_irt:irt_peptides_per_bin_nonlinear", 50, "Peptides sampled per bin for nonlinear iRT (0 = skip nonlinear)");
    defaults_.setMinInt("auto_irt:irt_peptides_per_bin_nonlinear", 0);
    // defaults_.setMaxInt("auto_irt:irt_peptides_per_bin_nonlinear", 10000);
    
    defaults_.setValue("auto_irt:linear_top_fraction", 0.4, "Top fraction of intense peptides to sample for linear iRT");
    defaults_.setMinFloat("auto_irt:linear_top_fraction", 0.01);
    defaults_.setMaxFloat("auto_irt:linear_top_fraction", 1.0);
    
    defaults_.setValue("auto_irt:nonlinear_top_fraction", 0.7, "Top fraction of intense peptides to sample for nonlinear iRT");
    defaults_.setMinFloat("auto_irt:nonlinear_top_fraction", 0.01);
    defaults_.setMaxFloat("auto_irt:nonlinear_top_fraction", 1.0);
    
    defaults_.setValue("auto_irt:irt_nonlinear_rt_extraction_window", 600.0, "Only extract RT around this value for non linear iRT calibration (set to -1 to use whole range)");
    defaults_.setMinFloat("auto_irt:irt_nonlinear_rt_extraction_window", -1.0);
    
    // Static iRT file parameters
    defaults_.setValue("files:linear_irt_file", "", "Path(s) to linear iRT transition file(s) (TraML, TSV, or PQP). Accepts a string of a single file path or multiple file paths (space-separated, 'run1_irt.pqp run2_irt.pqp ... runN_irt.pqp') for run-specific mapping (positional: nth entry -> nth run).");
    defaults_.setValue("files:nonlinear_irt_file", "", "Path(s) to nonlinear iRT transition file(s) (TraML, TSV, or PQP). Accepts a string of a single file path or multiple file paths (space-separated, 'run1_irt.pqp run2_irt.pqp ... runN_irt.pqp') for run-specific mapping (positional: nth entry -> nth run). Entries may be empty to indicate no nonlinear iRT for that run.");
    
    // Linear calibration parameters
    defaults_.setValue("linear:outlier_detection", "iter_residual", "Which outlier detection method to use for linear calibration (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
    defaults_.setValidStrings("linear:outlier_detection", {"iter_residual", "iter_jackknife", "ransac", "none"});
      
    // Nonlinear calibration parameters
    defaults_.setValue("nonlinear:outlier_detection", "iter_residual", "Which outlier detection method to use for nonlinear calibration (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
    defaults_.setValidStrings("nonlinear:outlier_detection", {"iter_residual", "iter_jackknife", "ransac", "none"});
    
    // Window estimation parameters
    defaults_.setValue("windows:estimate_rt", "true", "Estimate RT extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_rt", {"true", "false"});
    
    defaults_.setValue("windows:estimate_mz", "true", "Estimate m/z extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_mz", {"true", "false"});
    
    defaults_.setValue("windows:estimate_im", "true", "Estimate ion mobility extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_im", {"true", "false"});
    
    defaults_.setValue("windows:rt_percentile", 95.0, "Percentile for RT window estimation (25.0-99.9)");
    defaults_.setMinFloat("windows:rt_percentile", 25.0);
    defaults_.setMaxFloat("windows:rt_percentile", 99.9);
    
    // Window padding factors
    defaults_.setValue("windows:rt_estimation_padding_factor", 1.3, "A padding factor to multiply the estimated RT window by. For example, a factor of 1.3 will add a 30% padding to the estimated RT window, so if the estimated RT window is 144, then 43 will be added for a total estimated RT window of 187 seconds. A factor of 1.0 will not add any padding to the estimated window.");
    defaults_.setMinFloat("windows:rt_estimation_padding_factor", 1.0);
    
    // Quality control parameters
    defaults_.setValue("qc:min_rsq", 0.95, "Minimum R-squared required for RT peptides regression");
    defaults_.setMinFloat("qc:min_rsq", 0.0);
    defaults_.setMaxFloat("qc:min_rsq", 1.0);
    defaults_.setValue("qc:min_coverage", 0.6, "Minimum relative amount of RT peptides to keep");
    defaults_.setMinFloat("qc:min_coverage", 0.0);
    defaults_.setMaxFloat("qc:min_coverage", 1.0);

    // write defaults into Param object param_
    defaultsToParam_();
  }
  
  CalibrationWorkflow::~CalibrationWorkflow() = default;

  CalibrationWorkflow::CalibrationResult 
  CalibrationWorkflow::performCalibration(
    std::vector<OpenSwath::SwathMap>& swath_maps,
    OpenSwath::LightTargetedExperiment& transition_exp,
    ChromExtractParams& cp,
    ChromExtractParams& cp_ms1,
    const IrtExperiments& irt_experiments,
    const Param& feature_finder_param,
    const ChromExtractParams& cp_irt,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const Param& mrm_mapping_param,
    bool pasef,
    bool load_into_memory,
    const String& irt_trafo_out,
    const String& irt_mzml_out,
    Size debug_level)
  {
    // Validate that iRT experiments are ready
    if (!irt_experiments.is_prepared && irt_experiments.strategy != IrtStrategy::NULL_TRANSFORMATION)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "iRT experiments must be prepared before calling performCalibration. Call prepareIrtExperiments() first.");
    }
    
    // If the prepared experiments indicate there is a NULL_TRANSFORMATION strategy, skip calibration
    // and return a default/null calibration result (i.e., no RT transform will be applied).
    if (irt_experiments.strategy == IrtStrategy::NULL_TRANSFORMATION)
    {
      OPENMS_LOG_DEBUG << "IrtStrategy::NULL_TRANSFORMATION: skipping calibration and returning null transformation." << std::endl;
      return CalibrationResult();
    }

    CalibrationResult result;

    // Choose calibration workflow based on availability of nonlinear iRT data
    bool has_nonlinear_irt = !irt_experiments.nonlinear_irt.getTransitions().empty();

    if (has_nonlinear_irt)
    {
      result = performLinearThenNonlinearCalibration_(swath_maps, irt_experiments,
                                                     feature_finder_param, cp_irt, irt_detection_param,
                                                     calibration_param, mrm_mapping_param, pasef,
                                                     load_into_memory, irt_trafo_out, irt_mzml_out, debug_level);
    }
    else
    {
      result = performLinearCalibration_(swath_maps, irt_experiments,
                                        feature_finder_param, cp_irt, irt_detection_param,
                                        calibration_param, mrm_mapping_param, pasef,
                                        load_into_memory, irt_trafo_out, irt_mzml_out, debug_level);
    }
    
    // Apply IM correction back to the target library if needed (applies to both linear and nonlinear workflows)
    if (!result.im_trafo.getDataPoints().empty())
    {
      OPENMS_LOG_DEBUG << "Applying ion mobility correction to target library..." << std::endl;
      TransformationDescription im_trafo_inv = result.im_trafo;
      im_trafo_inv.invert();
      for (auto& compound : transition_exp.getCompounds())
      {
        if (compound.drift_time > 0)
        {
          compound.drift_time = im_trafo_inv.apply(compound.drift_time);
        }
      }
    }
    
    // Apply estimated extraction windows if requested
    applyEstimatedWindows_(result, cp, cp_ms1, pasef, 
                          feature_finder_param.getValue("use_ms1_ion_mobility").toBool());

    return result;
  }

  CalibrationWorkflow::CalibrationResult
  CalibrationWorkflow::performLinearCalibration_(
    std::vector<OpenSwath::SwathMap>& swath_maps,
    const IrtExperiments& irt_experiments,
    const Param& feature_finder_param,
    const ChromExtractParams& cp_irt,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const Param& mrm_mapping_param,
    bool pasef,
    bool load_into_memory,
    const String& irt_trafo_out,
    const String& irt_mzml_out,
    Size debug_level)
  {
    // Validate that we have linear iRT transitions (required for any calibration)
    if (irt_experiments.linear_irt.getTransitions().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "No linear iRT transitions available for calibration. CalibrationWorkflow requires at least linear iRT data to perform calibration.");
    }

    this->startProgress(0, 1, "Linear Calibration");

    CalibrationResult result;
      
    // Setup linear parameters with our outlier detection method
    Param linear_params = irt_detection_param;
    linear_params.setValue("outlierMethod", linear_outlier_detection_);

    // Configure calibration parameters - let SwathMapMassCorrection use its own m/z and IM parameters
    const Param& calibration_params_configured = calibration_param;

    TransformationDescription im_trafo;
    // Call the member implementation (moved into CalibrationWorkflow)
    result.rt_trafo = this->performRTNormalization(
      irt_experiments.linear_irt,
      swath_maps,
      im_trafo,
      min_rsq_,
      min_coverage_,
      feature_finder_param,
      cp_irt,
      linear_params,
      calibration_params_configured,
      mrm_mapping_param,
      irt_mzml_out,
      debug_level,
      pasef,
      load_into_memory);

    // Store the ion mobility transformation
    result.im_trafo = im_trafo;

    // Retrieve estimated windows (set by doDataNormalization_ during performRTNormalization())
    result.ms2_mz_window_ppm = this->estimated_mz_window_;
    result.ms2_im_window = this->estimated_im_window_;
    result.ms1_mz_window_ppm = this->estimated_ms1_mz_window_;
    result.ms1_im_window = this->estimated_ms1_im_window_;
    
    // Estimate RT window from transformation using configured parameters
    result.estimated_rt_window = result.rt_trafo.estimateWindow(
      windows_rt_percentile_ / 100.0,  // Convert percentage to fraction
      true,                           // Invert for RT units
      true,                           // Full width
      rt_estimation_padding_factor_); // User-configured padding

    // Save transformation if requested
    if (!irt_trafo_out.empty())
    {
      FileHandler().storeTransformations(irt_trafo_out, result.rt_trafo, 
                                        {FileTypes::TRANSFORMATIONXML});
    }
    this->endProgress();
    return result;
  }

  CalibrationWorkflow::CalibrationResult
  CalibrationWorkflow::performLinearThenNonlinearCalibration_(
    std::vector<OpenSwath::SwathMap>& swath_maps,
    const IrtExperiments& irt_experiments,
    const Param& feature_finder_param,
    const ChromExtractParams& cp_irt,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const Param& mrm_mapping_param,
    bool pasef,
    bool load_into_memory,
    const String& irt_trafo_out,
    const String& irt_mzml_out,
    Size debug_level)
  {
    // This method focuses on linear + nonlinear iRT-based calibration only
    
    // Step 1: Perform linear calibration (no m/z calibration to avoid double-application)
    
    Param linear_calibration_param = calibration_param;
    linear_calibration_param.setValue("mz_correction_function", "none");
    Param linear_detection_param = irt_detection_param;
    linear_detection_param.setValue("alignmentMethod", "linear");
    linear_detection_param.setValue("outlierMethod", linear_outlier_detection_);
    
    CalibrationResult linear_result = performLinearCalibration_(
      swath_maps, irt_experiments,
      feature_finder_param, cp_irt, linear_detection_param,
      linear_calibration_param, mrm_mapping_param, pasef,
      load_into_memory, irt_trafo_out, irt_mzml_out, debug_level);  
    
    // Step 2: Perform nonlinear refinement
    this->startProgress(0, 1, "Nonlinear Calibration");
    
    // Extract chromatograms for nonlinear iRT peptides using linear transformation
    std::vector<OpenMS::MSChromatogram> nl_chromatograms;
    ChromExtractParams cp_irt_nl = cp_irt;
    // Use potentially limited RT window for nonlinear extraction
    double nl_rt_window = (double)param_.getValue("auto_irt:irt_nonlinear_rt_extraction_window");
    if (nl_rt_window > 0)
    {
      cp_irt_nl.rt_extraction_window = nl_rt_window;
    }

    // Collect chromatograms using the default handler (same logic as previous helper)
    {
      std::unique_ptr<IChromatogramHandler> provider = IChromatogramHandler::createDefault();
      nl_chromatograms = provider->collectIrtChromatogramsForIrt(swath_maps, irt_experiments.nonlinear_irt, mrm_mapping_param, cp_irt_nl, linear_result.rt_trafo, pasef, load_into_memory);
    }
    
    // Setup nonlinear parameters
    Param nl_params = irt_detection_param;
    nl_params.setValue("estimateBestPeptides", "true"); // Enable outlier detection for nonlinear
    nl_params.setValue("outlierMethod", nonlinear_outlier_detection_); // Use nonlinear-specific outlier method
    
    // Configure calibration parameters - let SwathMapMassCorrection use its own m/z and IM parameters
    const Param& calibration_params_configured = calibration_param;
    
    TransformationDescription im_trafo;
    TransformationDescription nonlinear_trafo = this->doDataNormalization_(
      irt_experiments.nonlinear_irt,
      nl_chromatograms,
      im_trafo,
      swath_maps,
      min_rsq_,
      min_coverage_,
      feature_finder_param,
      nl_params,
      calibration_params_configured,
      pasef);

    // Prepare final result
    CalibrationResult final_result;
    final_result.rt_trafo = nonlinear_trafo;
    final_result.im_trafo = im_trafo;  // Store the ion mobility transformation
    final_result.ms2_mz_window_ppm = this->estimated_mz_window_;
    final_result.ms2_im_window = this->estimated_im_window_;
    final_result.ms1_mz_window_ppm = this->estimated_ms1_mz_window_;
    final_result.ms1_im_window = this->estimated_ms1_im_window_;
    
    final_result.estimated_rt_window = final_result.rt_trafo.estimateWindow(
      windows_rt_percentile_ / 100.0,  // Convert percentage to fraction
      true,                           // Invert for RT units  
      true,                           // Full width
      rt_estimation_padding_factor_); // User-configured padding

    // Save nonlinear transformation if requested
    if (!irt_trafo_out.empty())
    {
      String nonlinear_path = irt_trafo_out;
      const String ext = ".trafoXML";
      if (nonlinear_path.hasSuffix(ext))
      {
        nonlinear_path = nonlinear_path.substr(0, nonlinear_path.size() - ext.size());
        nonlinear_path += "_nonlinear.trafoXML";
      }
      FileHandler().storeTransformations(nonlinear_path, final_result.rt_trafo, 
                                       {FileTypes::TRANSFORMATIONXML});
    }


    this->endProgress();

    return final_result;
  }

  TransformationDescription CalibrationWorkflow::performRTNormalization(
    const OpenSwath::LightTargetedExperiment& irt_transitions,
    std::vector< OpenSwath::SwathMap > & swath_maps,
    TransformationDescription& im_trafo,
    double min_rsq,
    double min_coverage,
    const Param& feature_finder_param,
    const ChromExtractParams& cp_irt,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const Param& mrm_mapping_param,
    const String& irt_mzml_out,
    Size debug_level,
    bool pasef,
    bool load_into_memory)
  {
    std::vector< OpenMS::MSChromatogram > irt_chromatograms;
    TransformationDescription trafo; // dummy

    // collect & map chromatograms for iRT calibration.
    // The provider delegates to MRMChromHandler or DIAChromHandler as needed;
    {
      std::unique_ptr<IChromatogramHandler> provider = IChromatogramHandler::createDefault();
      irt_chromatograms = provider->collectIrtChromatogramsForIrt(swath_maps, irt_transitions, mrm_mapping_param, cp_irt, TransformationDescription(), pasef, load_into_memory);
    }

    // Summarize chromatogram extraction quality — empty chromatograms are a leading
    // indicator that downstream calibration will fail (insufficient RT bin coverage).
    // The per-batch "Detected N empty chromatograms" warnings from DIAChromHandler
    // scroll past easily; this one-line summary is the aggregated view.
    {
      Size nr_empty = std::count_if(irt_chromatograms.begin(), irt_chromatograms.end(),
        [](const MSChromatogram& c) { return c.empty(); });
      Size nr_total = irt_chromatograms.size();
      Size nr_nonempty = nr_total - nr_empty;
      if (nr_total > 0 && nr_empty > nr_total / 2)
      {
        OPENMS_LOG_WARN << "[OpenSwath] iRT chromatogram extraction: " << nr_nonempty << "/"
                        << nr_total << " non-empty (" << std::fixed << std::setprecision(1)
                        << (100.0 * nr_nonempty / nr_total) << "%). "
                        << "Over half are empty - iRT calibration may fail. "
                        << "Common causes: wrong SWATH boundaries, poor data quality, "
                        << "or unsorted m/z peaks in the input spectra." << std::endl;
      }
      else if (nr_total > 0)
      {
        OPENMS_LOG_INFO << "[OpenSwath] iRT chromatogram extraction: " << nr_nonempty << "/"
                        << nr_total << " non-empty (" << std::fixed << std::setprecision(1)
                        << (100.0 * nr_nonempty / nr_total) << "%)" << std::endl;
      }
    }

    // debug output of the iRT chromatograms
    String irt_mzml_out_local = irt_mzml_out;
    if (irt_mzml_out_local.empty() && debug_level > 1)
    {
      irt_mzml_out_local = "debug_irts.mzML";
    }
    if (!irt_mzml_out_local.empty())
    {
      try
      {
        PeakMap exp;
        exp.setChromatograms(irt_chromatograms);
        FileHandler().storeExperiment(irt_mzml_out_local, exp, {FileTypes::MZML});
      }
      catch (OpenMS::Exception::UnableToCreateFile& /*e*/)
      {
        OPENMS_LOG_DEBUG << "Error creating file " + irt_mzml_out_local + ", not writing out iRT chromatogram file"  << '\n';
      }
      catch (OpenMS::Exception::BaseException& /*e*/)
      {
        OPENMS_LOG_DEBUG << "Error writing to file " + irt_mzml_out_local + ", not writing out iRT chromatogram file"  << '\n';
      }
    }
    OPENMS_LOG_DEBUG << "Extracted number of chromatograms from iRT files: " << irt_chromatograms.size() <<  std::endl;

    // After collecting and optionally mapping iRT chromatograms, run the
    // data-normalization routine which performs peak picking and computes
    // the RT transformation. Return that transformation to the caller.
    TransformationDescription trafo_out = doDataNormalization_(irt_transitions,
                                                              irt_chromatograms,
                                                              im_trafo,
                                                              swath_maps,
                                                              min_rsq,
                                                              min_coverage,
                                                              feature_finder_param,
                                                              irt_detection_param,
                                                              calibration_param,
                                                              pasef);
    return trafo_out;
  }

  TransformationDescription CalibrationWorkflow::doDataNormalization_(
    const OpenSwath::LightTargetedExperiment& targeted_exp,
    const std::vector< OpenMS::MSChromatogram >& chromatograms,
    TransformationDescription& im_trafo,
    std::vector< OpenSwath::SwathMap > & swath_maps,
    double min_rsq,
    double min_coverage,
    const Param& default_ffparam,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const bool pasef)
  {
    bool estimateBestPeptides = irt_detection_param.getValue("estimateBestPeptides").toBool();

    // 1. Estimate the retention time range of the iRT peptides over all assays
    std::pair<double,double> RTRange = OpenSwathHelper::estimateRTRange(targeted_exp);

    // 2. Store the peptide retention times in an intermediate map
    std::unordered_map<OpenMS::String, double> PeptideRTMap;
    PeptideRTMap.reserve(targeted_exp.getCompounds().size());
    for (Size i = 0; i < targeted_exp.getCompounds().size(); i++)
    {
      PeptideRTMap[targeted_exp.getCompounds()[i].id] = targeted_exp.getCompounds()[i].rt;
    }

    // 3. Pick input chromatograms to identify RT pairs from the input data
    const OpenSwath::LightTargetedExperiment& transition_exp_used = targeted_exp;

    // Change the feature finding parameters:
    //  - no RT score (since we don't know the correct retention time)
    //  - no RT window
    //  - no elution model score
    //  - no peak quality (use all peaks)
    //  - if best peptides should be used, use peak quality
    MRMFeatureFinderScoring featureFinder;
    Param feature_finder_param(default_ffparam);
    feature_finder_param.setValue("Scores:use_rt_score", "false");
    feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    feature_finder_param.setValue("rt_extraction_window", -1.0);
    feature_finder_param.setValue("stop_report_after_feature", 1);
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 1.0); // set to 1.0 in all cases
    feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false"); // no peak quality -> take all peaks!

    double irt_mz_w = calibration_param.getValue("mz_extraction_window");
    bool irt_ppm = calibration_param.getValue("mz_extraction_window_ppm").toBool();
    feature_finder_param.setValue("irt_mz_extraction_window", irt_mz_w);
    feature_finder_param.setValue("irt_mz_extraction_window_unit", irt_ppm ? "ppm" : "Th");

    if (estimateBestPeptides)
    {
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", irt_detection_param.getValue("InitialQualityCutoff"));
    }
    featureFinder.setParameters(feature_finder_param);

    FeatureMap featureFile; // for results
    OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map; // for results
    std::vector<OpenSwath::SwathMap> empty_swath_maps;
    TransformationDescription empty_trafo; // empty transformation

    // Prepare the data with the chromatograms
    std::shared_ptr<PeakMap > xic_map(new PeakMap);
    xic_map->setChromatograms(chromatograms);
    OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(xic_map));

    featureFinder.setStrictFlag(false); // TODO remove this, it should be strict (e.g. all transitions need to be present for RT norm)
    featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp_used, empty_trafo, empty_swath_maps, transition_group_map);

    // 4. Find most likely correct feature for each compound and add it to the
    // "pairs" vector by computing pairs of iRT and real RT.
    //
    // Note that the quality threshold will only be applied if
    // estimateBestPeptides is true
    std::vector<std::pair<double, double> > pairs; // store the RT pairs to write the output trafoXML
    std::map<std::string, double> best_features = OpenSwathHelper::simpleFindBestFeature(transition_group_map,
      estimateBestPeptides, irt_detection_param.getValue("OverallQualityCutoff"));

    // Create pairs vector and store peaks
    std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> trgrmap_allpeaks; // store all peaks above cutoff
    for (std::map<std::string, double>::iterator it = best_features.begin(); it != best_features.end(); ++it)
    {
      pairs.emplace_back(it->second, PeptideRTMap[it->first]); // pair<exp_rt, theor_rt>
      auto tg_it = transition_group_map.find(it->first);
      if (tg_it != transition_group_map.end())
      {
        trgrmap_allpeaks[it->first] = &tg_it->second;
      }
    }

    // 5. Perform the outlier detection
    std::vector<std::pair<double, double> > pairs_corrected;
    String outlier_method = irt_detection_param.getValue("outlierMethod").toString();
    if (outlier_method == "iter_residual" || outlier_method == "iter_jackknife")
    {
      pairs_corrected = MRMRTNormalizer::removeOutliersIterative(pairs, min_rsq, min_coverage,
      irt_detection_param.getValue("useIterativeChauvenet").toBool(), outlier_method);
    }
    else if (outlier_method == "ransac")
    {
      // First, estimate of the maximum deviation from RT that is tolerated:
      //   Because 120 min gradient can have around 4 min elution shift, we use
      //   a default value of 3 % of the gradient to find upper RT threshold (3.6 min).
      double pcnt_rt_threshold = irt_detection_param.getValue("RANSACMaxPercentRTThreshold");
      double max_rt_threshold = (RTRange.second - RTRange.first) * pcnt_rt_threshold / 100.0;

      pairs_corrected = MRMRTNormalizer::removeOutliersRANSAC(pairs, min_rsq, min_coverage,
        irt_detection_param.getValue("RANSACMaxIterations"), max_rt_threshold,
        irt_detection_param.getValue("RANSACSamplingSize"));
    }
    else if (outlier_method == "none")
    {
      pairs_corrected = pairs;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        String("Illegal argument '") + outlier_method +
        "' used for outlierMethod (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none').");
    }

    // 6. Check whether the found peptides fulfill the binned coverage criteria
    // set by the user.
    if (estimateBestPeptides)
    {
      bool enoughPeptides = MRMRTNormalizer::computeBinnedCoverage(RTRange, pairs_corrected,
        irt_detection_param.getValue("NrRTBins"),
        irt_detection_param.getValue("MinPeptidesPerBin"),
        irt_detection_param.getValue("MinBinsFilled") );

      if (!enoughPeptides)
      {
        const int nr_bins = irt_detection_param.getValue("NrRTBins");
        const int min_pep = irt_detection_param.getValue("MinPeptidesPerBin");
        const int min_filled = irt_detection_param.getValue("MinBinsFilled");
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("iRT calibration failed: insufficient RT coverage after outlier removal.\n") +
          "  Chromatograms extracted:    " + String(chromatograms.size()) + "\n" +
          "  iRT peptides matched:       " + String(pairs.size()) + " (before outlier removal)\n" +
          "  After outlier removal:      " + String(pairs_corrected.size()) + " (exp_RT, theor_RT) pairs\n" +
          "  Binning requires:           " + String(min_filled) + "/" + String(nr_bins) + " RT bins "
          "with >= " + String(min_pep) + " peptides each\n" +
          "  RT range:                   " + String(RTRange.first, 1) + " - " + String(RTRange.second, 1) + " s\n" +
          "Common causes:\n"
          "  - Too few iRT/CiRT peptides in the transition library for the LC gradient\n"
          "  - Low data quality (check the 'Detected N empty chromatograms' warnings above)\n"
          "  - Wrong SWATH window boundaries (precursor m/z mismatch between library and data)");
      }
    }
    if (pairs_corrected.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        String("Less than 2 iRT normalization peptides after outlier removal (have ") +
        String(pairs_corrected.size()) + " from " + String(pairs.size()) + " initial matches, " +
        String(chromatograms.size()) + " chromatograms extracted). "
        "Not enough for an RT correction.");
    }

    // 7. Select the "correct" peaks for m/z (and IM) correction (e.g. remove those not
    // part of the linear regression)
    std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> trgrmap_final; // store all peaks above cutoff
    for (const auto& it : trgrmap_allpeaks)
    {
      if (it.second->getFeatures().empty() ) {continue;}
      const MRMFeature& feat = it.second->getBestFeature();

      // Check if the current feature is in the list of pairs used for the
      // linear RT regression (using other features may result in wrong
      // calibration values).
      // Matching only by RT is not perfect but should work for most cases.
      for (Size pit = 0; pit < pairs_corrected.size(); pit++)
      {
        if (fabs(feat.getRT() - pairs_corrected[pit].first ) < 1e-2)
        {
          trgrmap_final[ it.first ] = it.second;
          break;
        }
      }
    }

    // 8. Correct m/z (and IM) deviations using SwathMapMassCorrection
    // m/z correction is done with the -irt_im_extraction parameters
    SwathMapMassCorrection mc;
    mc.setParameters(calibration_param);

    mc.correctMZ(trgrmap_final, targeted_exp, swath_maps, pasef);
    mc.correctIM(trgrmap_final, targeted_exp, swath_maps, pasef, im_trafo);

    // Get estimated extraction windows (store in this workflow object)
    this->estimated_mz_window_ = mc.getFragmentMzWindow();
    this->estimated_im_window_ = mc.getFragmentImWindow();
    this->estimated_ms1_mz_window_ = mc.getPrecursorMzWindow();
    this->estimated_ms1_im_window_ = mc.getPrecursorImWindow();

    // 9. store RT transformation, using the selected model
    TransformationDescription trafo_out;
    trafo_out.setDataPoints(pairs_corrected);
    Param model_params;
    model_params.setValue("symmetric_regression", "false");
    model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
    model_params.setValue("auto_span", irt_detection_param.getValue("lowess:auto_span"));
    model_params.setValue("auto_span_min", irt_detection_param.getValue("lowess:auto_span_min"));
    model_params.setValue("auto_span_max", irt_detection_param.getValue("lowess:auto_span_max"));
    model_params.setValue("auto_span_grid", irt_detection_param.getValue("lowess:auto_span_grid"));
    model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
    String model_type = irt_detection_param.getValue("alignmentMethod").toString();
    trafo_out.fitModel(model_type, model_params);

    return trafo_out;
  }

  void CalibrationWorkflow::applyEstimatedWindows_(
    const CalibrationResult& result,
    ChromExtractParams& cp,
    ChromExtractParams& cp_ms1,
    bool pasef,
    bool use_ms1_im) const
  {
    // Apply RT window
    if (windows_estimate_rt_ && result.estimated_rt_window > 0)
    {
      applyWindow_("RT", result.estimated_rt_window, cp.rt_extraction_window, 
                  cp.rt_extraction_window, true, true);
      if (result.estimated_rt_window > 1000)
      {
        OPENMS_LOG_WARN << "Estimated RT extraction window is fairly large (" 
                        << result.estimated_rt_window 
                        << " seconds). If you are certain this is okay for your data, then ignore this warning. Otherwise, please verify that the calibration was successful by outputting the debugging calibration files or adjust the `windows:rt_percentile` to a lower value to restrict the residual distribution that the window is estimated from." 
                        << std::endl;
      }
    }
    
    // Apply MS2 m/z window (only if user is using ppm)
    if (windows_estimate_mz_ && result.ms2_mz_window_ppm > 0 && cp.ppm)
    {
      applyWindow_("MS2 m/z (ppm)", result.ms2_mz_window_ppm, cp.mz_extraction_window,
                  cp.mz_extraction_window, true, true);
    }
    
    // Apply MS2 ion mobility window
    if (windows_estimate_im_ && result.ms2_im_window > 0)
    {
      applyWindow_("MS2 ion mobility (1/k0)", result.ms2_im_window, cp.im_extraction_window,
                  cp.im_extraction_window, pasef, true);
    }
    
    // Apply MS1 m/z window (only if user is using ppm)  
    if (windows_estimate_mz_ && result.ms1_mz_window_ppm > 0 && cp_ms1.ppm)
    {
      applyWindow_("MS1 m/z (ppm)", result.ms1_mz_window_ppm, cp_ms1.mz_extraction_window,
                  cp_ms1.mz_extraction_window, true, true);
    }
    
    // Apply MS1 ion mobility window
    if (windows_estimate_im_ && result.ms1_im_window > 0)
    {
      applyWindow_("MS1 ion mobility (1/k0)", result.ms1_im_window, cp_ms1.im_extraction_window,
                  cp_ms1.im_extraction_window, pasef && use_ms1_im, true);
    }
  }

  void CalibrationWorkflow::applyWindow_(
    const char* label,
    double estimate,
    double& dst_param,
    double user_value,
    bool applicable,
    bool commit) const
  {
    if (!applicable)
    {
      OPENMS_LOG_INFO << "[Estimated] " << label
                      << " window: not applicable; keeping user value "
                      << user_value << std::endl;
      return;
    }

    if (!isValidWindow_(estimate))
    {
      OPENMS_LOG_WARN << "[Estimated] " << label
                      << " window estimate invalid (estimated=" << estimate
                      << "); keeping user value " << user_value << std::endl;
      return;
    }

    if (commit)
    {
      OPENMS_LOG_INFO << "[Estimated] " << label
                      << " window applied: " << estimate
                      << " (was " << user_value << ")" << std::endl;
      dst_param = estimate;
    }
    else
    {
      OPENMS_LOG_INFO << "[Estimated] " << label
                      << " window estimated: " << estimate
                      << "; keeping user value " << user_value << std::endl;
    }
  }

  bool CalibrationWorkflow::isValidWindow_(double v, double min_positive) const noexcept
  {
    return std::isfinite(v) && (v > min_positive);
  }



  OpenSwath::LightTargetedExperiment
  CalibrationWorkflow::loadIrtExperimentFromFile_(
    const String& irt_file_path,
    const String& label) const
  {
    OpenSwath::LightTargetedExperiment irt_exp;
    
    if (irt_file_path.empty())
    {
      return irt_exp; // Empty experiment if no file provided
    }
    
    if (!File::exists(irt_file_path))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, irt_file_path);
    }
    
    OPENMS_LOG_INFO << "Loading " << label << " iRT experiment from: " << irt_file_path << std::endl;
    
    try
    {
      FileTypes::Type file_type = FileHandler::getType(irt_file_path);
      
      if (file_type == FileTypes::TSV)
      {
        TransitionTSVFile tsv_reader;
        // Use default parameters for TSV reader
        TargetedExperiment targeted_exp;
        tsv_reader.convertTSVToTargetedExperiment(irt_file_path.c_str(), file_type, targeted_exp);
        OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, irt_exp);
      }
      else if (file_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_reader;
        TargetedExperiment targeted_exp;
        pqp_reader.convertPQPToTargetedExperiment(irt_file_path.c_str(), targeted_exp, false);
        OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, irt_exp);
      }
      else if (file_type == FileTypes::TRAML)
      {
        TraMLFile traml_reader;
        TargetedExperiment targeted_exp;
        traml_reader.load(irt_file_path, targeted_exp);
        OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, irt_exp);
      }
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        String("Unsupported iRT file format: ") + FileTypes::typeToName(file_type));
      }
      
      return irt_exp;
    }
    catch (const Exception::BaseException& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                String("Failed to load ") + label + " iRT file " + irt_file_path + ": " + e.getMessage(),
                                irt_file_path);
    }
  }


  IrtStrategy
  CalibrationWorkflow::determineIrtStrategy(
    const OpenSwath::LightTargetedExperiment& full_transition_exp,
    size_t num_runs) const
  {
    // Priority 1: If static iRT files are configured
    // If user supplied run-specific lists, prefer RUN_SPECIFIC when lists contain multiple entries
    if (!linear_irt_files_list_.empty())
    {
      if (linear_irt_files_list_.size() == 1)
      {
        OPENMS_LOG_DEBUG << "Single entry in files:linear_irt_file - treating as STATIC_FILES" << std::endl;
        return IrtStrategy::STATIC_FILES;
      }
      // multiple entries -> require num_runs match
      if (linear_irt_files_list_.size() != num_runs)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "files:linear_irt_file contains different number of entries (" + String((Int)linear_irt_files_list_.size()) + ") than the number of runs (" + String((Int)num_runs) + ") - provide one file per run or use STATIC_FILES");
      }
      // Validate nonlinear list length if present (must be same length or empty or size==1)
      if (!nonlinear_irt_files_list_.empty() && nonlinear_irt_files_list_.size() != linear_irt_files_list_.size() && nonlinear_irt_files_list_.size() != 1)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "files:nonlinear_irt_file length must be 0, 1, or equal to files:linear_irt_file length when using RUN_SPECIFIC");
      }
      OPENMS_LOG_DEBUG << "Using RUN_SPECIFIC iRT files (positional mapping) for " << num_runs << " runs." << std::endl;
      return IrtStrategy::RUN_SPECIFIC;
    }

    if (!linear_irt_file_.empty())
    {
      if (num_runs > 1)
      {
        OPENMS_LOG_DEBUG << "Static iRT files configured for " << num_runs << " runs - using STATIC_FILES strategy" << std::endl;
        return IrtStrategy::STATIC_FILES;
      }
      else
      {
        OPENMS_LOG_DEBUG << "Static iRT files configured for single run - using STATIC_FILES strategy" << std::endl;
        return IrtStrategy::STATIC_FILES;
      }
    }
    
    // Priority 2: If auto-iRT is enabled and we have transition library
    if (auto_irt_enabled_ && !full_transition_exp.getTransitions().empty())
    {
      if (num_runs > 1)
      {
        OPENMS_LOG_DEBUG << "Full transition library available for " << num_runs << " runs - using SAMPLE_ONCE strategy for consistency" << std::endl;
        return IrtStrategy::SAMPLE_ONCE;
      }
      else
      {
        OPENMS_LOG_DEBUG << "Full transition library available for single run - using SAMPLE_PER_RUN strategy" << std::endl;
        return IrtStrategy::SAMPLE_PER_RUN;
      }
    }
    
    // No iRT data available - return NONE strategy so callers can decide to skip calibration
    OPENMS_LOG_DEBUG << "No iRT data available: determineIrtStrategy() -> IrtStrategy::NULL_TRANSFORMATION" << std::endl;
    return IrtStrategy::NULL_TRANSFORMATION;
    }
  
  CalibrationWorkflow::IrtExperiments 
  CalibrationWorkflow::prepareIrtExperiments(
    IrtStrategy strategy,
    const OpenSwath::LightTargetedExperiment& full_transition_exp,
    const std::vector<String>& priority_peptides,
    size_t run_index,
    const IrtExperiments* cached_irts)
  {
    IrtExperiments result;
    result.strategy = strategy;
    result.is_prepared = false;
    
    switch (strategy) 
    {
      case IrtStrategy::STATIC_FILES:
        OPENMS_LOG_DEBUG << "Preparing static iRT experiments from configured files" << std::endl;
        
        // Load linear iRT experiment (required)
        if (!linear_irt_file_.empty())
        {
          result.linear_irt = loadIrtExperimentFromFile_(linear_irt_file_, "linear");
        }
        else
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "STATIC_FILES strategy requires linear_irt_file parameter to be configured");
        }
        
        // Load nonlinear iRT experiment (optional)
        if (!nonlinear_irt_file_.empty())
        {
          result.nonlinear_irt = loadIrtExperimentFromFile_(nonlinear_irt_file_, "nonlinear");
        }
        
        result.is_prepared = true;
        break;
        
      case IrtStrategy::SAMPLE_ONCE:
        if (cached_irts && cached_irts->is_prepared)
        {
          OPENMS_LOG_DEBUG << "Reusing cached sampled iRT experiments for run " << run_index << std::endl;
          result = *cached_irts; // Copy the cached experiments
        }
        else
        {
          OPENMS_LOG_DEBUG << "Sampling iRT experiments once from transition library (run " << run_index << ")" << std::endl;
          
          if (full_transition_exp.getTransitions().empty())
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "SAMPLE_ONCE strategy requires full_transition_exp to be provided");
          }
          
          // Convert vector<String> to unordered_set<string> for priority peptides
          std::unordered_set<std::string> priority_set;
          for (const auto& pep : priority_peptides) {
            priority_set.insert(pep.c_str());
          }
          
          // Generate linear iRT experiment
          result.linear_irt = OpenSwathHelper::sampleExperiment(
            full_transition_exp,
            auto_irt_irt_bins_,
            auto_irt_irt_peptides_per_bin_,
            auto_irt_irt_seed_,
            false,  // sort_by_intensity
            auto_irt_linear_top_fraction_,
            priority_set);
          
          // Generate nonlinear iRT experiment if requested
          if (auto_irt_irt_peptides_per_bin_nonlinear_ > 0)
          {
            result.nonlinear_irt = OpenSwathHelper::sampleExperiment(
              full_transition_exp,
              auto_irt_irt_bins_nonlinear_,
              auto_irt_irt_peptides_per_bin_nonlinear_,
              auto_irt_irt_seed_,
              false,  // sort_by_intensity
              auto_irt_nonlinear_top_fraction_,
              priority_set);
          }
          
          result.is_prepared = true;
        }
        break;
        
      case IrtStrategy::SAMPLE_PER_RUN:
      {
        OPENMS_LOG_DEBUG << "Sampling fresh iRT experiments for run " << run_index << std::endl;
        
        if (full_transition_exp.getTransitions().empty())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "SAMPLE_PER_RUN strategy requires full_transition_exp to be provided");
        }
        
        // Convert vector<String> to unordered_set<string> for priority peptides
        std::unordered_set<std::string> priority_set;
        for (const auto& pep : priority_peptides) {
          priority_set.insert(pep.c_str());
        }
        
        // Use run_index as additional seed variation for per-run sampling
        Size run_specific_seed = auto_irt_irt_seed_ + run_index;
        
        // Generate linear iRT experiment
        result.linear_irt = OpenSwathHelper::sampleExperiment(
          full_transition_exp,
          auto_irt_irt_bins_,
          auto_irt_irt_peptides_per_bin_,
          run_specific_seed,
          true,  // sort_by_intensity
          auto_irt_linear_top_fraction_,
          priority_set);
        
        // Generate nonlinear iRT experiment if requested
        if (auto_irt_irt_peptides_per_bin_nonlinear_ > 0)
        {
          result.nonlinear_irt = OpenSwathHelper::sampleExperiment(
            full_transition_exp,
            auto_irt_irt_bins_nonlinear_,
            auto_irt_irt_peptides_per_bin_nonlinear_,
            run_specific_seed,
            true,  // sort_by_intensity
            auto_irt_nonlinear_top_fraction_,
            priority_set);
        }
        
        result.is_prepared = true;
        break;
      }
        
      case IrtStrategy::RUN_SPECIFIC:
      {
        OPENMS_LOG_DEBUG << "Loading run-specific iRT experiments for run " << run_index << std::endl;
        // Positional mapping: nth entry in files lists corresponds to nth run (0-based)
        if (linear_irt_files_list_.empty())
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "RUN_SPECIFIC selected but files:linear_irt_file is empty");
        }
        if (run_index >= linear_irt_files_list_.size())
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Run index out of range for files:linear_irt_file (run_index=" + String((Int)run_index) + ", list_size=" + String((Int)linear_irt_files_list_.size()) + ")");
        }
        result.linear_irt = loadIrtExperimentFromFile_(linear_irt_files_list_[run_index], "linear");

        // Nonlinear file may be empty, single entry (apply to all) or positional
        if (!nonlinear_irt_files_list_.empty())
        {
          if (nonlinear_irt_files_list_.size() == 1)
          {
            // Single file apply to all runs
            if (!nonlinear_irt_files_list_[0].empty())
            {
              result.nonlinear_irt = loadIrtExperimentFromFile_(nonlinear_irt_files_list_[0], "nonlinear");
            }
          }
          else
          {
            if (run_index >= nonlinear_irt_files_list_.size())
            {
              throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "Run index out of range for files:nonlinear_irt_file (run_index=" + String((Int)run_index) + ", list_size=" + String((Int)nonlinear_irt_files_list_.size()) + ")");
            }
            if (!nonlinear_irt_files_list_[run_index].empty())
            {
              result.nonlinear_irt = loadIrtExperimentFromFile_(nonlinear_irt_files_list_[run_index], "nonlinear");
            }
          }
        }

        result.is_prepared = true;
        break;
      }

      case IrtStrategy::NULL_TRANSFORMATION:
      {
        OPENMS_LOG_DEBUG << "No iRT strategy selected - skipping preparation of iRT experiments" << std::endl;
        result.is_prepared = false;
        break;
      }
        
      default:
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Unknown IRT strategy provided");
    }
    
    return result;
  }

  void CalibrationWorkflow::updateMembers_()
  {
    // Static iRT file parameters
    linear_irt_file_ = param_.getValue("files:linear_irt_file").toString();
    nonlinear_irt_file_ = param_.getValue("files:nonlinear_irt_file").toString();
    
    // iRT peptide sampling parameters
    auto_irt_enabled_ = param_.getValue("auto_irt:enabled").toBool();
    auto_irt_irt_bins_ = (int)param_.getValue("auto_irt:irt_bins");
    auto_irt_irt_peptides_per_bin_ = (int)param_.getValue("auto_irt:irt_peptides_per_bin");
    auto_irt_irt_seed_ = (int)param_.getValue("auto_irt:irt_seed");
    auto_irt_irt_bins_nonlinear_ = (int)param_.getValue("auto_irt:irt_bins_nonlinear");
    auto_irt_irt_peptides_per_bin_nonlinear_ = (int)param_.getValue("auto_irt:irt_peptides_per_bin_nonlinear");
    auto_irt_linear_top_fraction_ = (double)param_.getValue("auto_irt:linear_top_fraction");
    auto_irt_nonlinear_top_fraction_ = (double)param_.getValue("auto_irt:nonlinear_top_fraction");
    
    // Linear calibration parameters
    linear_outlier_detection_ = param_.getValue("linear:outlier_detection").toString();
    
    // Nonlinear calibration parameters
    nonlinear_outlier_detection_ = param_.getValue("nonlinear:outlier_detection").toString();
    
    // Window estimation parameters
    windows_estimate_rt_ = param_.getValue("windows:estimate_rt").toBool();
    windows_estimate_mz_ = param_.getValue("windows:estimate_mz").toBool();
    windows_estimate_im_ = param_.getValue("windows:estimate_im").toBool();
    windows_rt_percentile_ = (double)param_.getValue("windows:rt_percentile");
    rt_estimation_padding_factor_ = (double)param_.getValue("windows:rt_estimation_padding_factor");
    
    // Quality control parameters
    min_rsq_ = (double)param_.getValue("qc:min_rsq");
    min_coverage_ = (double)param_.getValue("qc:min_coverage");

    // Run-specific iRT file lists
    linear_irt_files_list_.clear();
    nonlinear_irt_files_list_.clear();

    {
      StringList tmp;
      String legacy = param_.getValue("files:linear_irt_file").toString();
      legacy.split(' ', tmp);
      linear_irt_files_list_ = tmp;
    }

    {
      StringList tmp;
      String legacy = param_.getValue("files:nonlinear_irt_file").toString();
      legacy.split(' ', tmp);
      nonlinear_irt_files_list_ = tmp;
    }
  }

} // namespace OpenMS