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
    defaults_.setValue("auto_irt:enabled", "true", "Whether to sample iRTs on‐the‐fly (true) from the input targeted transition file (instead of passing specific iRT files). This may be useful if standard iRTs (Biognosys iRT kit) were not spiked-in. If set to false, and no additional iRT files are provided, and no transformation is provided, then no calibration is performed.");
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
    
    // === Static iRT file parameters ===
    defaults_.setValue("files:linear_irt_file", "", "Path to linear iRT transition file (TraML, TSV, or PQP)");
    defaults_.setValue("files:nonlinear_irt_file", "", "Path to nonlinear iRT transition file (TraML, TSV, or PQP)");
    
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
    defaults_.setValue("windows:estimate_rt", "true", "Estimate RT extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_rt", {"true", "false"});
    
    defaults_.setValue("windows:estimate_mz", "true", "Estimate m/z extraction windows from calibration");
    defaults_.setValidStrings("windows:estimate_mz", {"true", "false"});
    
    defaults_.setValue("windows:estimate_im", "true", "Estimate ion mobility extraction windows from calibration");
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
    
    // Window padding factors
    defaults_.setValue("windows:rt_estimation_padding_factor", 1.3, "A padding factor to multiply the estimated RT window by. For example, a factor of 1.3 will add a 30% padding to the estimated RT window, so if the estimated RT window is 144, then 43 will be added for a total estimated RT window of 187 seconds. A factor of 1.0 will not add any padding to the estimated window.");
    defaults_.setMinFloat("windows:rt_estimation_padding_factor", 1.0);
    
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
    const String& irt_mzml_out)
  {
    // Validate that iRT experiments are ready
    if (!irt_experiments.is_prepared)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "iRT experiments must be prepared before calling performCalibration. Call prepareIrtExperiments() first.");
    }

    CalibrationResult result;

    // Choose calibration workflow based on availability of nonlinear iRT data
    bool has_nonlinear_irt = !irt_experiments.nonlinear_irt.getTransitions().empty();

    if (has_nonlinear_irt)
    {
      OPENMS_LOG_INFO << "Performing linear + nonlinear calibration workflow..." << std::endl;
      result = performLinearThenNonlinearCalibration_(swath_maps, transition_exp, irt_experiments,
                                                     feature_finder_param, cp_irt, irt_detection_param,
                                                     calibration_param, mrm_mapping_param, pasef,
                                                     load_into_memory, irt_trafo_out, irt_mzml_out);
    }
    else
    {
      OPENMS_LOG_INFO << "Performing linear calibration workflow..." << std::endl;
      result = performLinearCalibration_(swath_maps, transition_exp, irt_experiments,
                                        feature_finder_param, cp_irt, irt_detection_param,
                                        calibration_param, mrm_mapping_param, pasef,
                                        load_into_memory, irt_trafo_out, irt_mzml_out);
    }
    
    // Apply estimated extraction windows if requested
    applyEstimatedWindows_(result, cp, cp_ms1, pasef, 
                          feature_finder_param.getValue("use_ms1_ion_mobility").toBool());

    OPENMS_LOG_INFO << "Calibration completed. Used " << result.num_irt_peptides_used 
                    << " iRT peptides (R²=" << result.rt_rsq << ", coverage=" 
                    << result.coverage_fraction << ")" << std::endl;

    return result;
  }

  CalibrationWorkflow::CalibrationResult
  CalibrationWorkflow::performLinearCalibration_(
    std::vector<OpenSwath::SwathMap>& swath_maps,
    OpenSwath::LightTargetedExperiment& transition_exp,
    const IrtExperiments& irt_experiments,
    const Param& feature_finder_param,
    const ChromExtractParams& cp_irt,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const Param& mrm_mapping_param,
    bool pasef,
    bool load_into_memory,
    const String& irt_trafo_out,
    const String& irt_mzml_out)
  {
    CalibrationResult result;
    
    // Note: This method focuses on iRT-based calibration only
    
    if (!irt_experiments.linear_irt.getTransitions().empty())
    {
      // Perform iRT-based calibration
      OPENMS_LOG_INFO << "Performing iRT-based linear calibration with " 
                      << irt_experiments.linear_irt.getTransitions().size() 
                      << " transitions..." << std::endl;
      
      // Create calibration workflow
      OpenSwathCalibrationWorkflow calibration_wf;
      calibration_wf.setLogType(getLogType());
      
      TransformationDescription im_trafo;
      result.rt_trafo = calibration_wf.performRTNormalization(
        irt_experiments.linear_irt, 
        swath_maps, 
        im_trafo,
        linear_min_rsq_, 
        min_coverage_,
        feature_finder_param,
        cp_irt,
        irt_detection_param,
        calibration_param,
        mrm_mapping_param,
        irt_mzml_out,
        0,  // debug_level
        pasef,
        load_into_memory);
        
      // Retrieve estimated windows
      result.ms2_mz_window_ppm = calibration_wf.getEstimatedMzWindow();
      result.ms2_im_window = calibration_wf.getEstimatedImWindow();
      result.ms1_mz_window_ppm = calibration_wf.getEstimatedMs1MzWindow();
      result.ms1_im_window = calibration_wf.getEstimatedMs1ImWindow();
      
      // Estimate RT window from transformation
      result.estimated_rt_window = result.rt_trafo.estimateWindow(
        0.99, true, true, 1.0);  // Use default padding factor

      // Save transformation if requested
      if (!irt_trafo_out.empty())
      {
        FileHandler().storeTransformations(irt_trafo_out, result.rt_trafo, 
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
    std::vector<OpenSwath::SwathMap>& swath_maps,
    OpenSwath::LightTargetedExperiment& transition_exp,
    const IrtExperiments& irt_experiments,
    const Param& feature_finder_param,
    const ChromExtractParams& cp_irt,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const Param& mrm_mapping_param,
    bool pasef,
    bool load_into_memory,
    const String& irt_trafo_out,
    const String& irt_mzml_out)
  {
    // This method focuses on linear + nonlinear iRT-based calibration only
    
    // Step 1: Perform linear calibration (no m/z calibration to avoid double-application)
    OPENMS_LOG_INFO << "Step 1: Performing linear calibration..." << std::endl;
    
    Param linear_calibration_param = calibration_param;
    linear_calibration_param.setValue("mz_correction_function", "none");
    Param linear_detection_param = irt_detection_param;
    linear_detection_param.setValue("alignmentMethod", "linear");
    
    CalibrationResult linear_result = performLinearCalibration_(
      swath_maps, transition_exp, irt_experiments,
      feature_finder_param, cp_irt, linear_detection_param,
      linear_calibration_param, mrm_mapping_param, pasef,
      load_into_memory, "", irt_mzml_out);  // Don't save intermediate trafo
    
    // Step 2: Perform nonlinear refinement
    OPENMS_LOG_INFO << "Step 2: Performing nonlinear calibration refinement..." << std::endl;
    
    OpenSwathCalibrationWorkflow nonlinear_wf;
    nonlinear_wf.setLogType(getLogType());
    
    // Extract chromatograms for nonlinear iRT peptides using linear transformation
    std::vector<OpenMS::MSChromatogram> nl_chromatograms;
    ChromExtractParams cp_irt_nl = cp_irt;
    // Use potentially limited RT window for nonlinear extraction
    // (this is often set in OpenSwathWorkflow via irt_nonlinear_rt_extraction_window)
    
    nonlinear_wf.simpleExtractChromatograms_(
      swath_maps,
      irt_experiments.nonlinear_irt,
      nl_chromatograms,
      linear_result.rt_trafo,  // Use linear result as starting point
      cp_irt_nl,
      mrm_mapping_param,
      pasef,
      load_into_memory);
    
    // Setup nonlinear parameters
    Param nl_params = irt_detection_param;
    nl_params.setValue("estimateBestPeptides", "true"); // Enable outlier detection
    
    TransformationDescription im_trafo;
    TransformationDescription nonlinear_trafo = nonlinear_wf.doDataNormalization_(
      irt_experiments.nonlinear_irt,
      nl_chromatograms,
      im_trafo,
      swath_maps,
      linear_min_rsq_,
      min_coverage_,
      feature_finder_param,
      nl_params,
      calibration_param,  // Use full calibration (with m/z correction)
      pasef);

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
      0.99, true, true, 1.0);  // Use default padding factor

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

    return final_result;
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


  IrtStrategy
  CalibrationWorkflow::determineIrtStrategy(
    const OpenSwath::LightTargetedExperiment& full_transition_exp,
    size_t num_runs) const
  {
    // Priority 1: If static iRT files are configured
    if (!linear_irt_file_.empty())
    {
      if (num_runs > 1)
      {
        OPENMS_LOG_INFO << "Static iRT files configured for " << num_runs << " runs - using STATIC_FILES strategy" << std::endl;
        return IrtStrategy::STATIC_FILES;
      }
      else
      {
        OPENMS_LOG_INFO << "Static iRT files configured for single run - using STATIC_FILES strategy" << std::endl;
        return IrtStrategy::STATIC_FILES;
      }
    }
    
    // Priority 2: If auto-iRT is enabled and we have transition library
    if (auto_irt_enabled_ && !full_transition_exp.getTransitions().empty())
    {
      if (num_runs > 1)
      {
        OPENMS_LOG_INFO << "Full transition library available for " << num_runs << " runs - using SAMPLE_ONCE strategy for consistency" << std::endl;
        return IrtStrategy::SAMPLE_ONCE;
      }
      else
      {
        OPENMS_LOG_INFO << "Full transition library available for single run - using SAMPLE_PER_RUN strategy" << std::endl;
        return IrtStrategy::SAMPLE_PER_RUN;
      }
    }
    
    // No iRT data available - CalibrationWorkflow requires iRT data
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "CalibrationWorkflow requires iRT data. Please provide either static iRT files (linear_irt_file parameter) or enable auto_irt and provide a full transition library for auto-sampling.");
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
        OPENMS_LOG_INFO << "Preparing static iRT experiments from configured files" << std::endl;
        
        // Load linear iRT experiment (required)
        if (!linear_irt_file_.empty())
        {
          result.linear_irt = loadIrtExperimentFromFile_(linear_irt_file_, "linear");
          OPENMS_LOG_INFO << "Loaded linear iRT experiment from file (" << result.linear_irt.getTransitions().size() << " transitions)" << std::endl;
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
          OPENMS_LOG_INFO << "Loaded nonlinear iRT experiment from file (" << result.nonlinear_irt.getTransitions().size() << " transitions)" << std::endl;
        }
        
        result.is_prepared = true;
        break;
        
      case IrtStrategy::SAMPLE_ONCE:
        if (cached_irts && cached_irts->is_prepared)
        {
          OPENMS_LOG_INFO << "Reusing cached sampled iRT experiments for run " << run_index << std::endl;
          result = *cached_irts; // Copy the cached experiments
        }
        else
        {
          OPENMS_LOG_INFO << "Sampling iRT experiments once from transition library (run " << run_index << ")" << std::endl;
          
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
          
          OPENMS_LOG_INFO << "Generated " << result.linear_irt.getTransitions().size() 
                          << " transitions for linear iRT calibration" << std::endl;
          
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
            
            OPENMS_LOG_INFO << "Generated " << result.nonlinear_irt.getTransitions().size()
                            << " transitions for nonlinear iRT calibration" << std::endl;
          }
          
          result.is_prepared = true;
        }
        break;
        
      case IrtStrategy::SAMPLE_PER_RUN:
      {
        OPENMS_LOG_INFO << "Sampling fresh iRT experiments for run " << run_index << std::endl;
        
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
          false,  // sort_by_intensity
          auto_irt_linear_top_fraction_,
          priority_set);
        
        OPENMS_LOG_INFO << "Generated " << result.linear_irt.getTransitions().size() 
                        << " transitions for linear iRT calibration" << std::endl;
        
        // Generate nonlinear iRT experiment if requested
        if (auto_irt_irt_peptides_per_bin_nonlinear_ > 0)
        {
          result.nonlinear_irt = OpenSwathHelper::sampleExperiment(
            full_transition_exp,
            auto_irt_irt_bins_nonlinear_,
            auto_irt_irt_peptides_per_bin_nonlinear_,
            run_specific_seed,
            false,  // sort_by_intensity
            auto_irt_nonlinear_top_fraction_,
            priority_set);
          
          OPENMS_LOG_INFO << "Generated " << result.nonlinear_irt.getTransitions().size()
                          << " transitions for nonlinear iRT calibration" << std::endl;
        }
        
        result.is_prepared = true;
        break;
      }
        
      case IrtStrategy::RUN_SPECIFIC:
      {
        OPENMS_LOG_INFO << "Loading run-specific iRT experiments for run " << run_index << std::endl;
        // TODO: Implement run-specific file loading logic
        // This would require run-specific file naming convention
        OPENMS_LOG_WARN << "RUN_SPECIFIC iRT preparation not yet implemented - using empty experiments" << std::endl;
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
    // === Static iRT file parameters ===
    linear_irt_file_ = param_.getValue("files:linear_irt_file").toString();
    nonlinear_irt_file_ = param_.getValue("files:nonlinear_irt_file").toString();
    
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
    
    // === Quality control parameters (from linear section for now) ===
    min_coverage_ = 0.6; // Default minimum coverage
    
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