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
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  CalibrationWorkflow::CalibrationWorkflow() :
    DefaultParamHandler("CalibrationWorkflow"),
    ProgressLogger()
  {
    // === iRT peptide sampling parameters ===
    defaults_.setValue("auto_irt:enabled", "false", "Enable automatic iRT peptide sampling for calibration");
    defaults_.setValidStrings("auto_irt:enabled", {"true", "false"});
    
    defaults_.setValue("auto_irt:irt_bins", 100, "Number of RT bins for linear iRT sampling");
    defaults_.setMinInt("auto_irt:irt_bins", 1);
    defaults_.setMaxInt("auto_irt:irt_bins", 10000);
    
    defaults_.setValue("auto_irt:irt_peptides_per_bin", 5, "Peptides sampled per bin for linear iRT");
    defaults_.setMinInt("auto_irt:irt_peptides_per_bin", 1);
    defaults_.setMaxInt("auto_irt:irt_peptides_per_bin", 1000);
    
    defaults_.setValue("auto_irt:irt_seed", 5489, "RNG seed for reproducible sampling (0 = non-deterministic)");
    defaults_.setMinInt("auto_irt:irt_seed", 0);
    
    defaults_.setValue("auto_irt:irt_bins_nonlinear", 2000, "Number of RT bins for nonlinear iRT sampling");
    defaults_.setMinInt("auto_irt:irt_bins_nonlinear", 1);
    defaults_.setMaxInt("auto_irt:irt_bins_nonlinear", 10000);
    
    defaults_.setValue("auto_irt:irt_peptides_per_bin_nonlinear", 50, "Peptides sampled per bin for nonlinear iRT (0 = skip nonlinear)");
    defaults_.setMinInt("auto_irt:irt_peptides_per_bin_nonlinear", 0);
    defaults_.setMaxInt("auto_irt:irt_peptides_per_bin_nonlinear", 1000);
    
    defaults_.setValue("auto_irt:linear_top_fraction", 0.4, "Top fraction of intense peptides to sample for linear iRT");
    defaults_.setMinFloat("auto_irt:linear_top_fraction", 0.01);
    defaults_.setMaxFloat("auto_irt:linear_top_fraction", 1.0);
    
    defaults_.setValue("auto_irt:nonlinear_top_fraction", 0.7, "Top fraction of intense peptides to sample for nonlinear iRT");
    defaults_.setMinFloat("auto_irt:nonlinear_top_fraction", 0.01);
    defaults_.setMaxFloat("auto_irt:nonlinear_top_fraction", 1.0);
    
    // === Linear calibration parameters ===
    defaults_.setValue("linear:enabled", "true", "Perform linear RT calibration");
    defaults_.setValidStrings("linear:enabled", {"true", "false"});
    
    defaults_.setValue("linear:outlier_detection", "iter_residual", "Outlier detection method");
    defaults_.setValidStrings("linear:outlier_detection", {"iter_residual", "iter_jackknife", "ransac", "none"});
    
    defaults_.setValue("linear:min_rsq", 0.95, "Minimum R-squared required for linear calibration");
    defaults_.setMinFloat("linear:min_rsq", 0.0);
    defaults_.setMaxFloat("linear:min_rsq", 1.0);
    
    // === Nonlinear calibration parameters ===
    defaults_.setValue("nonlinear:enabled", "false", "Perform nonlinear RT calibration after linear");
    defaults_.setValidStrings("nonlinear:enabled", {"true", "false"});
    
    defaults_.setValue("nonlinear:method", "lowess", "Nonlinear fitting method");
    defaults_.setValidStrings("nonlinear:method", {"lowess", "bspline", "polynomial"});
    
    defaults_.setValue("nonlinear:asymmetric", "false", "Use asymmetric regression (huber loss)");
    defaults_.setValidStrings("nonlinear:asymmetric", {"true", "false"});
    
    defaults_.setValue("nonlinear:span", 0.75, "Lowess span parameter (0.2-2.0)");
    defaults_.setMinFloat("nonlinear:span", 0.2);
    defaults_.setMaxFloat("nonlinear:span", 2.0);
    
    // === Window estimation parameters ===
    defaults_.setValue("windows:estimate_rt", "false", "Estimate RT extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_rt", {"true", "false"});
    
    defaults_.setValue("windows:estimate_mz", "false", "Estimate m/z extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_mz", {"true", "false"});
    
    defaults_.setValue("windows:estimate_im", "false", "Estimate ion mobility extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_im", {"true", "false"});
    
    defaults_.setValue("windows:rt_percentile", 95.0, "Percentile for RT window estimation (80.0-99.9)");
    defaults_.setMinFloat("windows:rt_percentile", 80.0);
    defaults_.setMaxFloat("windows:rt_percentile", 99.9);
    
    defaults_.setValue("windows:mz_percentile", 95.0, "Percentile for m/z window estimation (80.0-99.9)");
    defaults_.setMinFloat("windows:mz_percentile", 80.0);
    defaults_.setMaxFloat("windows:mz_percentile", 99.9);
    
    defaults_.setValue("windows:im_percentile", 95.0, "Percentile for IM window estimation (80.0-99.9)");
    defaults_.setMinFloat("windows:im_percentile", 80.0);
    defaults_.setMaxFloat("windows:im_percentile", 99.9);
    
    // === Quality control parameters ===
    defaults_.setValue("qc:fail_on_insufficient_peptides", "true", "Fail if insufficient peptides found");
    defaults_.setValidStrings("qc:fail_on_insufficient_peptides", {"true", "false"});
    
    defaults_.setValue("qc:fail_on_poor_fit", "true", "Fail if calibration fit quality is poor");
    defaults_.setValidStrings("qc:fail_on_poor_fit", {"true", "false"});
    
    defaults_.setValue("qc:fail_on_low_coverage", "true", "Fail if coverage fraction too low");
    defaults_.setValidStrings("qc:fail_on_low_coverage", {"true", "false"});

    // write defaults into Param object param_
    defaultsToParam_();
  }
  
  CalibrationWorkflow::~CalibrationWorkflow() = default;

  CalibrationWorkflow::CalibrationResult 
  CalibrationWorkflow::performCalibration(
    const String& trafo_in,
    std::vector<OpenSwath::SwathMap>& swath_maps,
    OpenSwath::LightTargetedExperiment& transition_exp,
    ChromExtractParams& cp,
    ChromExtractParams& cp_ms1,
    const CalibrationConfig& config)
  {
    // Validate configuration
    if (!validateConfig(config, trafo_in))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Invalid calibration configuration provided.");
    }

    CalibrationResult result;

    // Priority 1: Load RT transformation file if provided (rt_norm parameter)
    if (!trafo_in.empty())
    {
      OPENMS_LOG_INFO << "Loading user-provided RT transformation from: " << trafo_in << std::endl;
      try 
      {
        FileHandler().loadTransformations(trafo_in, result.rt_trafo, true);
        OPENMS_LOG_INFO << "Successfully loaded RT transformation with " 
                        << result.rt_trafo.getDataPoints().size() 
                        << " data points" << std::endl;
        
        // When user provides transformation, we don't estimate windows
        // (user is bypassing the iRT calibration entirely)
        result.ms2_mz_window_ppm = -1.0;  // Invalid/unused
        result.ms2_im_window = -1.0;      // Invalid/unused  
        result.ms1_mz_window_ppm = -1.0;  // Invalid/unused
        result.ms1_im_window = -1.0;      // Invalid/unused
        result.estimated_rt_window = -1.0; // Invalid/unused
        
        OPENMS_LOG_INFO << "Skipping iRT calibration and window estimation (user-provided transformation)" << std::endl;
        
        // Apply estimated extraction windows if requested (though estimates will be -1)
        applyEstimatedWindows_(result, cp, cp_ms1, config, config.pasef, 
                              config.feature_finder_param.getValue("use_ms1_ion_mobility").toBool());
        
        OPENMS_LOG_INFO << "RT transformation loaded successfully." << std::endl;
        return result;
      }
      catch (const Exception::BaseException& e)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  String("Failed to load RT transformation file ") + trafo_in + ": " + e.getMessage(),
                                  trafo_in);
      }
    }

    // Priority 2-4: Proceed with iRT-based calibration (file-based, pre-loaded, or auto-sampled)
    OPENMS_LOG_INFO << "No RT transformation file provided - proceeding with iRT-based calibration" << std::endl;

    // Set up the actual iRT experiments to use - prioritize files over pre-loaded experiments
    OpenSwath::LightTargetedExperiment linear_irt_exp = config.linear_irt_exp;
    OpenSwath::LightTargetedExperiment nonlinear_irt_exp = config.nonlinear_irt_exp;
    
    // Load from files if provided (takes precedence over pre-loaded experiments)
    if (!config.linear_irt_file.empty())
    {
      linear_irt_exp = loadIrtExperimentFromFile_(config.linear_irt_file, config.irt_tsv_reader_param, "linear");
    }
    if (!config.nonlinear_irt_file.empty())
    {
      nonlinear_irt_exp = loadIrtExperimentFromFile_(config.nonlinear_irt_file, config.irt_tsv_reader_param, "nonlinear");
    }
    
    // Auto-generate iRT experiments if enabled and full experiment provided
    if (param_.getValue("auto_irt:enabled").toBool() && 
        !config.full_transition_exp.getTransitions().empty())
    {
      OPENMS_LOG_INFO << "Auto-generating iRT experiments from " 
                      << config.full_transition_exp.getTransitions().size()
                      << " total transitions..." << std::endl;
      
      // Parameters for linear iRT sampling
      Size irt_bins = (Size)param_.getValue("auto_irt:irt_bins");
      Size irt_peptides_per_bin = (Size)param_.getValue("auto_irt:irt_peptides_per_bin");
      Size irt_seed = (Size)param_.getValue("auto_irt:irt_seed");
      double linear_top_fraction = (double)param_.getValue("auto_irt:linear_top_fraction");
      
      // Convert vector<String> to unordered_set<string> for priority peptides
      std::unordered_set<std::string> priority_set;
      for (const auto& pep : config.priority_peptides) {
        priority_set.insert(pep.c_str());
      }
      
      // Generate linear iRT experiment (or enhance existing one)
      if (linear_irt_exp.getTransitions().empty())
      {
        linear_irt_exp = OpenSwathHelper::sampleExperiment(
          config.full_transition_exp,
          irt_bins,
          irt_peptides_per_bin,
          irt_seed,
          false,  // sort_by_intensity
          linear_top_fraction,
          priority_set);
        
        OPENMS_LOG_INFO << "Generated " << linear_irt_exp.getTransitions().size() 
                        << " transitions for linear iRT calibration" << std::endl;
      }
      
      // Parameters for nonlinear iRT sampling
      Size irt_bins_nonlinear = (Size)param_.getValue("auto_irt:irt_bins_nonlinear");
      Size irt_peptides_per_bin_nonlinear = (Size)param_.getValue("auto_irt:irt_peptides_per_bin_nonlinear");
      double nonlinear_top_fraction = (double)param_.getValue("auto_irt:nonlinear_top_fraction");
      
      // Only generate nonlinear iRT if requested (peptides_per_bin > 0)
      if (irt_peptides_per_bin_nonlinear > 0 && nonlinear_irt_exp.getTransitions().empty())
      {
        nonlinear_irt_exp = OpenSwathHelper::sampleExperiment(
          config.full_transition_exp,
          irt_bins_nonlinear,
          irt_peptides_per_bin_nonlinear,
          irt_seed,
          false,  // sort_by_intensity
          nonlinear_top_fraction,
          priority_set);
        
        OPENMS_LOG_INFO << "Generated " << nonlinear_irt_exp.getTransitions().size()
                        << " transitions for nonlinear iRT calibration" << std::endl;
      }
    }

    // Create a working config with the final iRT experiments
    CalibrationConfig working_config = config;
    working_config.linear_irt_exp = linear_irt_exp;
    working_config.nonlinear_irt_exp = nonlinear_irt_exp;

    // Determine calibration strategy: linear-only vs linear+nonlinear
    bool has_nonlinear_irt = !working_config.nonlinear_irt_exp.getTransitions().empty();
    
    if (has_nonlinear_irt)
    {
      OPENMS_LOG_INFO << "Performing linear + nonlinear calibration workflow..." << std::endl;
      result = performLinearThenNonlinearCalibration_("", swath_maps, 
                                                     transition_exp, working_config);
    }
    else
    {
      OPENMS_LOG_INFO << "Performing linear calibration workflow..." << std::endl;
      result = performLinearCalibration_("", swath_maps, transition_exp, working_config);
    }
    
    // Apply estimated extraction windows if requested
    applyEstimatedWindows_(result, cp, cp_ms1, working_config, working_config.pasef, 
                          working_config.feature_finder_param.getValue("use_ms1_ion_mobility").toBool());

    OPENMS_LOG_INFO << "Calibration completed. Used " << result.num_irt_peptides_used 
                    << " iRT peptides (R²=" << result.rt_rsq << ", coverage=" 
                    << result.coverage_fraction << ")" << std::endl;

    return result;
  }

  CalibrationWorkflow::CalibrationResult
  CalibrationWorkflow::performLinearCalibration_(
    const String& trafo_in,
    std::vector<OpenSwath::SwathMap>& swath_maps,
    OpenSwath::LightTargetedExperiment& transition_exp,
    const CalibrationConfig& config)
  {
    CalibrationResult result;
    
    // Note: trafo_in is handled at the performCalibration level
    // This method focuses on iRT-based calibration only
    
    if (!config.linear_irt_exp.getTransitions().empty())
    {
      // Perform iRT-based calibration
      OPENMS_LOG_INFO << "Performing iRT-based linear calibration with " 
                      << config.linear_irt_exp.getTransitions().size() 
                      << " transitions..." << std::endl;
      
      // Create calibration workflow
      OpenSwathCalibrationWorkflow calibration_wf;
      calibration_wf.setLogType(getLogType());
      
      TransformationDescription im_trafo;
      result.rt_trafo = calibration_wf.performRTNormalization(
        config.linear_irt_exp, 
        swath_maps, 
        im_trafo,
        config.min_rsq, 
        config.min_coverage,
        config.feature_finder_param,
        config.cp_irt,
        config.irt_detection_param,
        config.calibration_param,
        config.mrm_mapping_param,
        config.irt_mzml_out,
        config.debug_level,
        config.pasef,
        config.load_into_memory);
        
      // Retrieve estimated windows
      result.ms2_mz_window_ppm = calibration_wf.getEstimatedMzWindow();
      result.ms2_im_window = calibration_wf.getEstimatedImWindow();
      result.ms1_mz_window_ppm = calibration_wf.getEstimatedMs1MzWindow();
      result.ms1_im_window = calibration_wf.getEstimatedMs1ImWindow();
      
      // Estimate RT window from transformation
      result.estimated_rt_window = result.rt_trafo.estimateWindow(
        0.99, true, true, config.rt_estimation_padding_factor);

      // Save transformation if requested
      if (!config.irt_trafo_out.empty())
      {
        FileHandler().storeTransformations(config.irt_trafo_out, result.rt_trafo, 
                                         {FileTypes::TRANSFORMATIONXML});
      }
    }
    else
    {
      OPENMS_LOG_INFO << "No iRT transitions provided, using identity transformation." << std::endl;
      // Create identity transformation
      // result.rt_trafo = TransformationDescription(); // identity
    }

    return result;
  }

  CalibrationWorkflow::CalibrationResult
  CalibrationWorkflow::performLinearThenNonlinearCalibration_(
    const String& trafo_in,
    std::vector<OpenSwath::SwathMap>& swath_maps,
    OpenSwath::LightTargetedExperiment& transition_exp,
    const CalibrationConfig& config)
  {
    // Note: trafo_in is handled at the performCalibration level
    // This method focuses on linear + nonlinear iRT-based calibration only
    
    // Step 1: Perform linear calibration (no m/z calibration to avoid double-application)
    OPENMS_LOG_INFO << "Step 1: Performing linear calibration..." << std::endl;
    
    CalibrationConfig linear_config = config;
    linear_config.calibration_param.setValue("mz_correction_function", "none");
    linear_config.irt_detection_param.setValue("alignmentMethod", "linear");
    
    CalibrationResult linear_result = performLinearCalibration_(
      "", swath_maps, transition_exp, linear_config);
    
    // Step 2: Perform nonlinear refinement
    OPENMS_LOG_INFO << "Step 2: Performing nonlinear calibration refinement..." << std::endl;
    
    OpenSwathCalibrationWorkflow nonlinear_wf;
    nonlinear_wf.setLogType(getLogType());
    
    // Extract chromatograms for nonlinear iRT peptides using linear transformation
    std::vector<OpenMS::MSChromatogram> nl_chromatograms;
    ChromExtractParams cp_irt_nl = config.cp_irt;
    // Use potentially limited RT window for nonlinear extraction
    // (this is often set in OpenSwathWorkflow via irt_nonlinear_rt_extraction_window)
    
    nonlinear_wf.simpleExtractChromatograms_(
      swath_maps,
      config.nonlinear_irt_exp,
      nl_chromatograms,
      linear_result.rt_trafo,  // Use linear result as starting point
      cp_irt_nl,
      config.mrm_mapping_param,
      config.pasef,
      config.load_into_memory);
    
    // Setup nonlinear parameters
    Param nl_params = config.irt_detection_param;
    nl_params.setValue("estimateBestPeptides", "true"); // Enable outlier detection
    
    TransformationDescription im_trafo;
    TransformationDescription nonlinear_trafo = nonlinear_wf.doDataNormalization_(
      config.nonlinear_irt_exp,
      nl_chromatograms,
      im_trafo,
      swath_maps,
      config.min_rsq,
      config.min_coverage,
      config.feature_finder_param,
      nl_params,
      config.calibration_param,  // Use full calibration (with m/z correction)
      config.pasef);

    // Apply IM correction back to the target library if needed
    if (!im_trafo.getDataPoints().empty())
    {
      TransformationDescription im_trafo_inv = im_trafo;
      im_trafo_inv.invert();
      for (auto& compound : transition_exp.getCompounds())
      {
        if (compound.drift_time > 0)
        {
          compound.drift_time = im_trafo_inv.apply(compound.drift_time);
        }
      }
    }
    
    // Prepare final result
    CalibrationResult final_result;
    final_result.rt_trafo = nonlinear_trafo;
    final_result.ms2_mz_window_ppm = nonlinear_wf.getEstimatedMzWindow();
    final_result.ms2_im_window = nonlinear_wf.getEstimatedImWindow();
    final_result.ms1_mz_window_ppm = nonlinear_wf.getEstimatedMs1MzWindow();
    final_result.ms1_im_window = nonlinear_wf.getEstimatedMs1ImWindow();
    
    final_result.estimated_rt_window = final_result.rt_trafo.estimateWindow(
      0.99, true, true, config.rt_estimation_padding_factor);

    // Save nonlinear transformation if requested
    if (!config.irt_trafo_out.empty())
    {
      String nonlinear_path = config.irt_trafo_out;
      const String ext = ".trafoXML";
      if (nonlinear_path.hasSuffix(ext))
      {
        nonlinear_path = nonlinear_path.substr(0, nonlinear_path.size() - ext.size());
        nonlinear_path += "_nonlinear.trafoXML";
      }
      FileHandler().storeTransformations(nonlinear_path, final_result.rt_trafo, 
                                       {FileTypes::TRANSFORMATIONXML});
    }

    return final_result;
  }

  void CalibrationWorkflow::applyEstimatedWindows_(
    const CalibrationResult& result,
    ChromExtractParams& cp,
    ChromExtractParams& cp_ms1,
    const CalibrationConfig& config,
    bool pasef,
    bool use_ms1_im) const
  {
    const auto& est_windows = config.use_estimated_windows;
    
    // Apply RT window
    if (est_windows.rt && result.estimated_rt_window > 0)
    {
      applyWindow_("RT", result.estimated_rt_window, cp.rt_extraction_window, 
                  cp.rt_extraction_window, true, true);
    }
    
    // Apply MS2 m/z window (only if user is using ppm)
    if (est_windows.mz && result.ms2_mz_window_ppm > 0 && cp.ppm)
    {
      applyWindow_("MS2 m/z (ppm)", result.ms2_mz_window_ppm, cp.mz_extraction_window,
                  cp.mz_extraction_window, true, true);
    }
    
    // Apply MS2 ion mobility window
    if (est_windows.im && result.ms2_im_window > 0)
    {
      applyWindow_("MS2 ion mobility (1/k0)", result.ms2_im_window, cp.im_extraction_window,
                  cp.im_extraction_window, pasef, true);
    }
    
    // Apply MS1 m/z window (only if user is using ppm)  
    if (est_windows.mz && result.ms1_mz_window_ppm > 0 && cp_ms1.ppm)
    {
      applyWindow_("MS1 m/z (ppm)", result.ms1_mz_window_ppm, cp_ms1.mz_extraction_window,
                  cp_ms1.mz_extraction_window, true, true);
    }
    
    // Apply MS1 ion mobility window
    if (est_windows.im && result.ms1_im_window > 0)
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

  bool CalibrationWorkflow::validateConfig(const CalibrationConfig& config) const
  {
    // Basic parameter validation
    if (config.min_rsq < 0.0 || config.min_rsq > 1.0)
    {
      OPENMS_LOG_ERROR << "Invalid min_rsq: " << config.min_rsq 
                       << " (must be between 0.0 and 1.0)" << std::endl;
      return false;
    }
    
    if (config.min_coverage < 0.0 || config.min_coverage > 1.0)
    {
      OPENMS_LOG_ERROR << "Invalid min_coverage: " << config.min_coverage 
                       << " (must be between 0.0 and 1.0)" << std::endl;
      return false;
    }
    
    // Check that we have some calibration method
    bool has_linear_irt = !config.linear_irt_exp.getTransitions().empty() || !config.linear_irt_file.empty();
    bool has_nonlinear_irt = !config.nonlinear_irt_exp.getTransitions().empty() || !config.nonlinear_irt_file.empty();
    bool has_full_exp_for_sampling = !config.full_transition_exp.getTransitions().empty();
    bool auto_irt_enabled = param_.getValue("auto_irt:enabled").toBool();
    
    if (!has_linear_irt && !has_nonlinear_irt && !(auto_irt_enabled && has_full_exp_for_sampling))
    {
      OPENMS_LOG_ERROR << "No calibration method available: provide iRT experiments or enable auto-iRT with full transition experiment" << std::endl;
      return false;
    }
    
    if (!has_linear_irt && has_nonlinear_irt && !(auto_irt_enabled && has_full_exp_for_sampling))
    {
      OPENMS_LOG_ERROR << "Cannot perform nonlinear calibration without linear iRT transitions" << std::endl;
      return false;
    }
    
    return true;
  }

  bool CalibrationWorkflow::validateConfig(const CalibrationConfig& config, const String& trafo_in) const
  {
    // Basic parameter validation (same as regular validateConfig)
    if (config.min_rsq < 0.0 || config.min_rsq > 1.0)
    {
      OPENMS_LOG_ERROR << "Invalid min_rsq: " << config.min_rsq 
                       << " (must be between 0.0 and 1.0)" << std::endl;
      return false;
    }
    
    if (config.min_coverage < 0.0 || config.min_coverage > 1.0)
    {
      OPENMS_LOG_ERROR << "Invalid min_coverage: " << config.min_coverage 
                       << " (must be between 0.0 and 1.0)" << std::endl;
      return false;
    }
    
    // Check that we have some calibration method - now including trafo_in
    bool has_trafo_file = !trafo_in.empty();
    bool has_linear_irt = !config.linear_irt_exp.getTransitions().empty() || !config.linear_irt_file.empty();
    bool has_nonlinear_irt = !config.nonlinear_irt_exp.getTransitions().empty() || !config.nonlinear_irt_file.empty();
    bool has_full_exp_for_sampling = !config.full_transition_exp.getTransitions().empty();
    bool auto_irt_enabled = param_.getValue("auto_irt:enabled").toBool();
    
    // Valid calibration methods:
    // 1. User-provided transformation file (rt_norm)
    // 2. iRT-based calibration (linear and/or nonlinear)  
    // 3. Auto-iRT sampling from full experiment
    if (!has_trafo_file && !has_linear_irt && !has_nonlinear_irt && !(auto_irt_enabled && has_full_exp_for_sampling))
    {
      OPENMS_LOG_ERROR << "No calibration method available: provide RT transformation file, iRT experiments, or enable auto-iRT with full transition experiment" << std::endl;
      return false;
    }
    
    // If user provides transformation file, validate it exists
    if (has_trafo_file)
    {
      if (!File::exists(trafo_in))
      {
        OPENMS_LOG_ERROR << "RT transformation file does not exist: " << trafo_in << std::endl;
        return false;
      }
      
      // When transformation file is provided, other calibration methods are ignored
      if (has_linear_irt || has_nonlinear_irt || (auto_irt_enabled && has_full_exp_for_sampling))
      {
        OPENMS_LOG_WARN << "RT transformation file provided - ignoring iRT calibration settings" << std::endl;
      }
    }
    else
    {
      // No transformation file - validate iRT-based calibration
      if (!has_linear_irt && has_nonlinear_irt && !(auto_irt_enabled && has_full_exp_for_sampling))
      {
        OPENMS_LOG_ERROR << "Cannot perform nonlinear calibration without linear iRT transitions" << std::endl;
        return false;
      }
    }
    
    return true;
  }

  OpenSwath::LightTargetedExperiment
  CalibrationWorkflow::loadIrtExperimentFromFile_(
    const String& irt_file_path,
    const Param& irt_tsv_reader_param,
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
        tsv_reader.setParameters(irt_tsv_reader_param);
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
      
      OPENMS_LOG_INFO << "Loaded " << irt_exp.getTransitions().size() << " transitions for " 
                      << label << " iRT calibration from " << irt_file_path << std::endl;
      
      return irt_exp;
    }
    catch (const Exception::BaseException& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                String("Failed to load ") + label + " iRT file " + irt_file_path + ": " + e.getMessage(),
                                irt_file_path);
    }
  }

  CalibrationWorkflow::CalibrationConfig 
  CalibrationWorkflow::getDefaultConfig()
  {
    CalibrationConfig config;
    
    // Quality control defaults
    config.min_rsq = 0.95;
    config.min_coverage = 0.6;
    
    // iRT file parameters defaults
    config.irt_tsv_reader_param = TransitionTSVFile().getDefaults();
    
    // Algorithm parameters (would typically be populated from TOPP tool parameters)
    config.feature_finder_param = MRMFeatureFinderScoring().getDefaults();
    config.cp_irt = ChromExtractParams(); // Default extraction parameters
    // ... other defaults would be set here
    
    // Settings defaults
    config.pasef = false;
    config.load_into_memory = false;
    config.debug_level = 0;
    
    // Window estimation defaults
    config.use_estimated_windows = {false, false, false}; // Conservative default
    config.rt_estimation_padding_factor = 1.3;
    config.im_estimation_padding_factor = 1.0;
    config.mz_estimation_padding_factor = 1.0;
    
    return config;
  }

  void CalibrationWorkflow::updateMembers_()
  {
    // === iRT peptide sampling parameters ===
    auto_irt_enabled_ = param_.getValue("auto_irt:enabled").toBool();
    auto_irt_irt_bins_ = (int)param_.getValue("auto_irt:irt_bins");
    auto_irt_irt_peptides_per_bin_ = (int)param_.getValue("auto_irt:irt_peptides_per_bin");
    auto_irt_irt_seed_ = (int)param_.getValue("auto_irt:irt_seed");
    auto_irt_irt_bins_nonlinear_ = (int)param_.getValue("auto_irt:irt_bins_nonlinear");
    auto_irt_irt_peptides_per_bin_nonlinear_ = (int)param_.getValue("auto_irt:irt_peptides_per_bin_nonlinear");
    auto_irt_linear_top_fraction_ = (double)param_.getValue("auto_irt:linear_top_fraction");
    auto_irt_nonlinear_top_fraction_ = (double)param_.getValue("auto_irt:nonlinear_top_fraction");
    
    // === Linear calibration parameters ===
    linear_enabled_ = param_.getValue("linear:enabled").toBool();
    linear_outlier_detection_ = param_.getValue("linear:outlier_detection").toString();
    linear_min_rsq_ = (double)param_.getValue("linear:min_rsq");
    
    // === Nonlinear calibration parameters ===
    nonlinear_enabled_ = param_.getValue("nonlinear:enabled").toBool();
    nonlinear_method_ = param_.getValue("nonlinear:method").toString();
    nonlinear_asymmetric_ = param_.getValue("nonlinear:asymmetric").toBool();
    nonlinear_span_ = (double)param_.getValue("nonlinear:span");
    
    // === Window estimation parameters ===
    windows_estimate_rt_ = param_.getValue("windows:estimate_rt").toBool();
    windows_estimate_mz_ = param_.getValue("windows:estimate_mz").toBool();
    windows_estimate_im_ = param_.getValue("windows:estimate_im").toBool();
    windows_rt_percentile_ = (double)param_.getValue("windows:rt_percentile");
    windows_mz_percentile_ = (double)param_.getValue("windows:mz_percentile");
    windows_im_percentile_ = (double)param_.getValue("windows:im_percentile");
    
    // === Quality control parameters ===
    qc_fail_on_insufficient_peptides_ = param_.getValue("qc:fail_on_insufficient_peptides").toBool();
    qc_fail_on_poor_fit_ = param_.getValue("qc:fail_on_poor_fit").toBool();
    qc_fail_on_low_coverage_ = param_.getValue("qc:fail_on_low_coverage").toBool();
  }

} // namespace OpenMS