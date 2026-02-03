// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/APPLICATIONS/OpenSwathBase.h>

namespace OpenMS
{
  /**
    @brief Orchestrates calibration workflows for OpenSWATH analysis.
    
    This class provides a unified interface for performing RT, m/z, and ion mobility
    calibration workflows. It supports both linear-only calibration and 
    linear+nonlinear calibration strategies, making calibration logic reusable
    across different OpenSWATH tools (DIA/SWATH, SRM/MRM, PRM).
    
    The orchestrator handles:
    - Parameter validation and setup
    - iRT peptide sampling (when auto_irt is enabled)
    - Linear calibration workflow
    - Nonlinear calibration workflow (optional)
    - Extraction window estimation and application
    - Result aggregation and reporting
    
    @htmlinclude OpenMS_CalibrationWorkflow.parameters
    
    @ingroup OpenSWATH
  */
  class OPENMS_DLLAPI CalibrationWorkflow : 
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    
    /// Selection flags for using auto-estimated extraction windows
    struct EstimateWindowsChoice
    {
      bool rt{false};   ///< Use estimated RT extraction window
      bool mz{false};   ///< Use estimated m/z extraction window  
      bool im{false};   ///< Use estimated ion mobility extraction window
    };
    
    /// Configuration for calibration workflow
    struct CalibrationConfig 
    {
      // === iRT Experiments ===
      /// iRT transitions for linear calibration (pre-loaded)
      OpenSwath::LightTargetedExperiment linear_irt_exp;
      /// iRT transitions for nonlinear calibration (pre-loaded, optional)
      OpenSwath::LightTargetedExperiment nonlinear_irt_exp;
      
      // === iRT File Paths (alternative to pre-loaded experiments) ===
      /// Path to linear iRT file (TraML/TSV/PQP format, empty = use linear_irt_exp)
      String linear_irt_file;
      /// Path to nonlinear iRT file (TraML/TSV/PQP format, empty = use nonlinear_irt_exp)
      String nonlinear_irt_file;
      /// Parameters for loading iRT files
      Param irt_tsv_reader_param;
      
      // === Auto-iRT Sampling (optional) ===
      /// Full transition experiment for auto-iRT sampling (empty = use provided iRT experiments)
      OpenSwath::LightTargetedExperiment full_transition_exp;
      /// Priority peptides for iRT sampling (empty = no priorities)
      std::vector<String> priority_peptides;
      
      // === Quality Control Parameters ===
      /// Minimum R-squared value for RT regression
      double min_rsq{0.95};
      /// Minimum coverage of chromatographic space
      double min_coverage{0.6};
      
      // === Algorithm Parameters ===
      /// Parameters for MRMFeatureFinderScoring
      Param feature_finder_param;
      /// Parameters for chromatogram extraction (iRT peptides)
      ChromExtractParams cp_irt;
      /// Parameters for iRT detection and outlier removal
      Param irt_detection_param;
      /// Parameters for m/z and IM calibration
      Param calibration_param;
      /// Parameters for MRM chromatogram mapping
      Param mrm_mapping_param;
      
      // === Acquisition Settings ===
      /// Whether data is PASEF (ion mobility) data
      bool pasef{false};
      /// Whether to load data into memory for processing
      bool load_into_memory{false};
      /// Debug level for output (0=minimal, higher=more verbose)
      Size debug_level{0};
      
      // === Output Files (optional) ===
      /// Output file for RT transformation (empty = no output)
      String irt_trafo_out;
      /// Output file for iRT chromatograms (empty = no output)  
      String irt_mzml_out;
      
      // === Window Estimation Settings ===
      /// Which extraction windows to auto-estimate and apply
      EstimateWindowsChoice use_estimated_windows;
      /// Padding factor for estimated RT windows (e.g., 1.3 = +30%)
      double rt_estimation_padding_factor{1.3};
      /// Padding factor for estimated IM windows
      double im_estimation_padding_factor{1.0};
      /// Padding factor for estimated m/z windows
      double mz_estimation_padding_factor{1.0};
    };
    
    /// Results from calibration workflow
    struct CalibrationResult 
    {
      // === Transformation ===
      /// RT normalization transformation (fitted)
      TransformationDescription rt_trafo;
      
      // === Estimated Extraction Windows ===
      /// MS2 m/z extraction window (full width, ppm). -1 if not computed.
      double ms2_mz_window_ppm{-1.0};
      /// MS2 ion mobility extraction window (full width). -1 if not applicable.
      double ms2_im_window{-1.0};
      /// MS1 m/z extraction window (full width, ppm). -1 if not computed.
      double ms1_mz_window_ppm{-1.0};
      /// MS1 ion mobility extraction window (full width). -1 if not applicable.
      double ms1_im_window{-1.0};
      /// Estimated RT extraction window (full width, seconds)
      double estimated_rt_window{-1.0};
      
      // === Quality Metrics ===
      /// Number of iRT peptides used in final calibration
      Size num_irt_peptides_used{0};
      /// R-squared of final RT transformation
      double rt_rsq{0.0};
      /// Coverage fraction achieved
      double coverage_fraction{0.0};
    };

    /// Default constructor
    CalibrationWorkflow();
    
    /// Destructor  
    ~CalibrationWorkflow() override;

    /**
      @brief Perform calibration workflow orchestration.
      
      This is the main entry point that orchestrates the entire calibration workflow.
      It automatically chooses between linear-only or linear+nonlinear calibration
      based on the configuration provided. This method handles both the simple
      linear calibration case and the complex linear+nonlinear case, exactly
      matching the logic flow in OpenSwathWorkflow.
      
      @param[in] trafo_in User-provided transformation file (empty = perform calibration)
      @param[in,out] swath_maps Raw SWATH/SRM data maps (modified in-place by calibration)
      @param[in,out] transition_exp Target transition experiment (IM values may be corrected)
      @param[in,out] cp Extraction parameters (windows may be updated with estimates)
      @param[in,out] cp_ms1 MS1 extraction parameters (windows may be updated)
      @param[in] config Calibration configuration and parameters
      
      @return Calibration results including transformations and estimated windows
      
      @throws Exception::IllegalArgument If configuration is invalid
      @throws Exception::MissingInformation If required iRT data is missing
    */
    CalibrationResult performCalibration(
      const String& trafo_in,
      std::vector<OpenSwath::SwathMap>& swath_maps,
      OpenSwath::LightTargetedExperiment& transition_exp,
      ChromExtractParams& cp,
      ChromExtractParams& cp_ms1,
      const CalibrationConfig& config);

    /**
      @brief Validate calibration configuration.
      
      Checks that the provided configuration is valid and all required
      parameters are present.
      
      @param[in] config Configuration to validate
      @return True if valid, false otherwise
    */
    bool validateConfig(const CalibrationConfig& config) const;
    
    /**
      @brief Validate configuration with optional transformation file.
      
      Extended validation that considers both the config and an optional
      RT transformation file. This allows validation when a user provides
      a pre-computed transformation file (rt_norm parameter).
      
      @param[in] config Configuration to validate  
      @param[in] trafo_in Path to RT transformation file (empty = not provided)
      @return True if valid, false otherwise
    */
    bool validateConfig(const CalibrationConfig& config, const String& trafo_in) const;
    
    /**
      @brief Load iRT transition experiment from file if path is provided.
      
      Helper method to load iRT experiments from file paths in CalibrationConfig.
      Supports TraML, TSV, and PQP formats.
      
      @param[in] irt_file_path Path to iRT file (empty = skip loading)
      @param[in] irt_tsv_reader_param Parameters for TSV reader
      @param[in] label Label for logging (e.g., "linear", "nonlinear")
      @return Loaded iRT experiment (empty if file path was empty)
      @throws Exception::FileNotFound if file doesn't exist
      @throws Exception::ParseError if file format is invalid
    */
    OpenSwath::LightTargetedExperiment loadIrtExperimentFromFile_(
      const String& irt_file_path,
      const Param& irt_tsv_reader_param,
      const String& label) const;
    
    /**
      @brief Get the default calibration configuration.
      
      Returns a default configuration that can be customized by the caller.
      
      @return Default calibration configuration
    */
    static CalibrationConfig getDefaultConfig();

  private:
    
    /// @name Parameter handling
    //@{
    /**
      @brief Update member variables from parameters.
      
      This method is called automatically when parameters change.
      Updates internal member variables from current parameter values.
    */
    void updateMembers_() override;
    //@}
    
    /// @name Private member variables for parameter caching  
    //@{
    bool auto_irt_enabled_;
    int auto_irt_irt_bins_;
    int auto_irt_irt_peptides_per_bin_;
    int auto_irt_irt_seed_;
    int auto_irt_irt_bins_nonlinear_;
    int auto_irt_irt_peptides_per_bin_nonlinear_;
    double auto_irt_linear_top_fraction_;
    double auto_irt_nonlinear_top_fraction_;
    
    bool linear_enabled_;
    String linear_outlier_detection_;
    double linear_min_rsq_;
    
    bool nonlinear_enabled_;
    String nonlinear_method_;
    bool nonlinear_asymmetric_;
    double nonlinear_span_;
    
    bool windows_estimate_rt_;
    bool windows_estimate_mz_;
    bool windows_estimate_im_;
    double windows_rt_percentile_;
    double windows_mz_percentile_;
    double windows_im_percentile_;
    
    bool qc_fail_on_insufficient_peptides_;
    bool qc_fail_on_poor_fit_;
    bool qc_fail_on_low_coverage_;
    //@}
    
    /// @name Private implementation methods
    //@{
    
    /**
      @brief Perform linear-only calibration workflow.
      
      @param[in] trafo_in User-provided transformation file
      @param[in,out] swath_maps SWATH data maps
      @param[in,out] transition_exp Target transitions  
      @param[in] config Calibration configuration
      
      @return Linear calibration results
    */
    CalibrationResult performLinearCalibration_(
      const String& trafo_in,
      std::vector<OpenSwath::SwathMap>& swath_maps,
      OpenSwath::LightTargetedExperiment& transition_exp,
      const CalibrationConfig& config);
    
    /**
      @brief Perform linear + nonlinear calibration workflow.
      
      First performs a linear calibration, then applies nonlinear refinement
      using a separate set of iRT transitions.
      
      @param[in] trafo_in User-provided transformation file  
      @param[in,out] swath_maps SWATH data maps
      @param[in,out] transition_exp Target transitions
      @param[in] config Calibration configuration
      
      @return Combined calibration results
    */
    CalibrationResult performLinearThenNonlinearCalibration_(
      const String& trafo_in,
      std::vector<OpenSwath::SwathMap>& swath_maps,
      OpenSwath::LightTargetedExperiment& transition_exp,
      const CalibrationConfig& config);
      
    /**
      @brief Apply estimated extraction windows to parameters.
      
      Updates the provided extraction parameters with auto-estimated windows
      based on calibration results and user preferences.
      
      @param[in] result Calibration results containing estimated windows
      @param[in,out] cp MS2 extraction parameters to update
      @param[in,out] cp_ms1 MS1 extraction parameters to update
      @param[in] config Configuration specifying which windows to apply
      @param[in] pasef Whether this is PASEF data (for IM applicability)
      @param[in] use_ms1_im Whether MS1 uses ion mobility
    */
    void applyEstimatedWindows_(
      const CalibrationResult& result,
      ChromExtractParams& cp,
      ChromExtractParams& cp_ms1, 
      const CalibrationConfig& config,
      bool pasef,
      bool use_ms1_im) const;
      
    /**
      @brief Validate and log an auto-estimated extraction window.
      
      @param[in] label Human-readable label for logging
      @param[in] estimate Estimated window value  
      @param[in,out] dst_param Parameter to update
      @param[in] user_value Current user-specified value
      @param[in] applicable Whether this window type is applicable
      @param[in] commit Whether to actually apply the estimate
    */
    void applyWindow_(
      const char* label,
      double estimate,
      double& dst_param,
      double user_value,
      bool applicable = true,
      bool commit = true) const;
      
    /**
      @brief Check if an estimated window value is valid.
      
      @param[in] v Window value to check
      @param[in] min_positive Minimum positive threshold
      @return True if window is valid for use
    */
    inline bool isValidWindow_(double v, double min_positive = 1e-9) const noexcept;
    //@}
  };

} // namespace OpenMS