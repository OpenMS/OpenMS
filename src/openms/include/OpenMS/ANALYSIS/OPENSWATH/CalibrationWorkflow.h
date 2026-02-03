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
    @brief Strategy for obtaining iRT (indexed Retention Time) calibrant peptides across multiple runs.
    
    This enum defines the different approaches for providing iRT peptides
    for retention time calibration, with consideration for multi-run consistency:
    
    - STATIC_FILES: Use pre-defined iRT peptide libraries from files (same for all runs). If all samples have spiked in iRTs, or common expected peptides.
    - SAMPLE_ONCE: Sample iRT peptides once from transition library, reuse for all runs  
    - SAMPLE_PER_RUN: Sample fresh iRT peptides independently for each run
    - RUN_SPECIFIC: Use different iRT files per run (each run has its own iRT files). These are usually generated from accompanying pseudo-DDA spectra analysis
  */
  enum class IrtStrategy 
  {
    STATIC_FILES,     
    SAMPLE_ONCE,        
    SAMPLE_PER_RUN,   
    RUN_SPECIFIC      
  };

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
    
    /// Prepared iRT experiments ready for calibration
    struct IrtExperiments
    {
      /// Linear iRT experiment (always required for calibration)
      OpenSwath::LightTargetedExperiment linear_irt;
      /// Nonlinear iRT experiment (optional, empty if not used)  
      OpenSwath::LightTargetedExperiment nonlinear_irt;
      /// Strategy used to prepare these experiments
      IrtStrategy strategy{IrtStrategy::STATIC_FILES};
      /// Whether experiments are available and ready for calibration
      bool is_prepared{false};
    };

    /// Default constructor
    CalibrationWorkflow();
    
    /// Destructor  
    ~CalibrationWorkflow() override;
    
    /**
      @brief Determine the appropriate IRT strategy based on available data.
      
      Analyzes the available iRT data to determine which IRT strategy
      should be used for calibration. Priority order:
      1. If static iRT files are configured → STATIC_FILES or RUN_SPECIFIC
      2. If full transition library is provided → SAMPLE_ONCE or SAMPLE_PER_RUN
      3. Otherwise → throws exception (CalibrationWorkflow requires iRT data)
      
      @param[in] full_transition_exp Full transition experiment for auto-sampling (empty if not available)
      @param[in] num_runs Total number of runs to process (affects strategy choice)
      @return The determined IRT strategy
      @throws Exception::MissingInformation if no iRT data is available
    */
    IrtStrategy determineIrtStrategy(
      const OpenSwath::LightTargetedExperiment& full_transition_exp,
      size_t num_runs = 1
    ) const;
    
    /**
      @brief Prepare iRT experiments based on the determined strategy.
      
      This function handles all aspects of iRT experiment preparation:
      - Loading static iRT files  
      - Auto-sampling from transition libraries with priority peptides
      - Validation and quality checks
      - Multi-run consistency management
      
      The prepared experiments are ready for use with performCalibration().
      
      @param[in] strategy The IRT strategy to use
      @param[in] full_transition_exp Full transition experiment for auto-sampling (required for sampling strategies)
      @param[in] priority_peptides Priority peptide sequences for sampling (empty = no priorities)
      @param[in] run_index Current run index (0-based, for run-specific strategies) 
      @param[in] cached_irts Previously prepared iRT experiments (for reusing sampled iRTs across runs)
      @return Prepared iRT experiments ready for calibration
      @throws Exception::MissingInformation if required data is not available
      @throws Exception::InvalidParameter if configuration is invalid
    */
    IrtExperiments prepareIrtExperiments(
      IrtStrategy strategy,
      const OpenSwath::LightTargetedExperiment& full_transition_exp,
      const std::vector<String>& priority_peptides,
      size_t run_index = 0,
      const IrtExperiments* cached_irts = nullptr
    );

    /**
      @brief Perform calibration workflow with pre-prepared iRT experiments.
      
      This is the main calibration entry point that accepts pre-prepared iRT experiments
      and performs RT, m/z, and ion mobility calibration. It automatically chooses 
      between linear-only or linear+nonlinear calibration based on whether nonlinear
      iRT experiments are provided.
      
      @param[in,out] swath_maps Raw SWATH/SRM data maps (modified in-place by calibration)
      @param[in,out] transition_exp Target transition experiment (IM values may be corrected)
      @param[in,out] cp Extraction parameters (windows may be updated with estimates)
      @param[in,out] cp_ms1 MS1 extraction parameters (windows may be updated)
      @param[in] irt_experiments Pre-prepared iRT experiments for calibration
      @param[in] feature_finder_param Parameters for MRMFeatureFinderScoring
      @param[in] cp_irt Extraction parameters for iRT peptides
      @param[in] irt_detection_param Parameters for iRT detection and outlier removal
      @param[in] calibration_param Parameters for m/z and IM calibration  
      @param[in] mrm_mapping_param Parameters for MRM chromatogram mapping
      @param[in] pasef Whether data is PASEF (ion mobility) data
      @param[in] load_into_memory Whether to load data into memory for processing
      @param[in] irt_trafo_out Output file for RT transformation (empty = no output)
      @param[in] irt_mzml_out Output file for iRT chromatograms (empty = no output)
      
      @return Calibration results including transformations and estimated windows
      
      @throws Exception::IllegalArgument If configuration is invalid
      @throws Exception::MissingInformation If iRT experiments are not ready
    */
    CalibrationResult performCalibration(
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
      bool pasef = false,
      bool load_into_memory = false,
      const String& irt_trafo_out = "",
      const String& irt_mzml_out = ""
    );

    /**
      @brief Load iRT transition experiment from file.
      
      Helper method to load iRT experiments from configured file paths.
      Supports TraML, TSV, and PQP formats.
      
      @param[in] irt_file_path Path to iRT file (empty = skip loading)
      @param[in] label Label for logging (e.g., "linear", "nonlinear")
      @return Loaded iRT experiment (empty if file path was empty)
      @throws Exception::FileNotFound if file doesn't exist
      @throws Exception::ParseError if file format is invalid
    */
    OpenSwath::LightTargetedExperiment loadIrtExperimentFromFile_(
      const String& irt_file_path,
      const String& label) const;

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
    
    // === iRT File Parameters ===
    String linear_irt_file_;                    ///< Path to linear iRT file
    String nonlinear_irt_file_;                 ///< Path to nonlinear iRT file
    
    // === Quality Control Parameters ===
    double min_rsq_;                            ///< Minimum R-squared value for RT regression
    double min_coverage_;                       ///< Minimum coverage of chromatographic space
    
    // === Auto-iRT Sampling Parameters ===
    bool auto_irt_enabled_;                     ///< Enable auto-iRT sampling
    int auto_irt_irt_bins_;                     ///< Number of RT bins for linear iRT sampling
    int auto_irt_irt_peptides_per_bin_;         ///< Peptides per bin for linear iRT
    int auto_irt_irt_seed_;                     ///< RNG seed for sampling
    int auto_irt_irt_bins_nonlinear_;           ///< Number of RT bins for nonlinear iRT
    int auto_irt_irt_peptides_per_bin_nonlinear_; ///< Peptides per bin for nonlinear iRT
    double auto_irt_linear_top_fraction_;       ///< Top fraction for linear sampling
    double auto_irt_nonlinear_top_fraction_;    ///< Top fraction for nonlinear sampling
    
    // === Linear Calibration Parameters ===
    bool linear_enabled_;                       ///< Enable linear calibration
    String linear_outlier_detection_;           ///< Outlier detection method
    double linear_min_rsq_;                     ///< Min R-squared for linear fit
    
    // === Nonlinear Calibration Parameters ===
    bool nonlinear_enabled_;                    ///< Enable nonlinear calibration
    String nonlinear_method_;                   ///< Nonlinear method name
    bool nonlinear_asymmetric_;                 ///< Use asymmetric nonlinear fit
    double nonlinear_span_;                     ///< Span parameter for nonlinear
    
    // === Window Estimation Parameters ===
    bool windows_estimate_rt_;                  ///< Estimate RT windows
    bool windows_estimate_mz_;                  ///< Estimate m/z windows
    bool windows_estimate_im_;                  ///< Estimate IM windows
    double windows_rt_percentile_;              ///< RT percentile for estimation
    double windows_mz_percentile_;              ///< m/z percentile for estimation
    double windows_im_percentile_;              ///< IM percentile for estimation
    double rt_estimation_padding_factor_;       ///< RT padding factor
    double im_estimation_padding_factor_;       ///< IM padding factor
    double mz_estimation_padding_factor_;       ///< m/z padding factor
    
    // === Quality Control Parameters ===
    bool qc_fail_on_insufficient_peptides_;     ///< Fail on insufficient peptides
    bool qc_fail_on_poor_fit_;                  ///< Fail on poor fit quality
    bool qc_fail_on_low_coverage_;              ///< Fail on low coverage
    
    //@}
    
    /// @name Private implementation methods
    //@{
    
    /**
      @brief Perform linear-only calibration workflow.
      
      @param[in,out] swath_maps SWATH data maps
      @param[in,out] transition_exp Target transitions  
      @param[in] irt_experiments Prepared iRT experiments
      @param[in] feature_finder_param Parameters for MRMFeatureFinderScoring
      @param[in] cp_irt Extraction parameters for iRT peptides
      @param[in] irt_detection_param Parameters for iRT detection
      @param[in] calibration_param Parameters for calibration
      @param[in] mrm_mapping_param Parameters for MRM mapping
      @param[in] pasef Whether this is PASEF data
      @param[in] load_into_memory Whether to load data into memory
      @param[in] irt_trafo_out Output transformation file
      @param[in] irt_mzml_out Output iRT chromatograms file
      
      @return Linear calibration results
    */
    CalibrationResult performLinearCalibration_(
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
      const String& irt_mzml_out);
    
    /**
      @brief Perform linear + nonlinear calibration workflow.
      
      First performs a linear calibration, then applies nonlinear refinement
      using a separate set of iRT transitions.
      
      @param[in,out] swath_maps SWATH data maps
      @param[in,out] transition_exp Target transitions
      @param[in] irt_experiments Prepared iRT experiments
      @param[in] feature_finder_param Parameters for MRMFeatureFinderScoring
      @param[in] cp_irt Extraction parameters for iRT peptides
      @param[in] irt_detection_param Parameters for iRT detection
      @param[in] calibration_param Parameters for calibration
      @param[in] mrm_mapping_param Parameters for MRM mapping
      @param[in] pasef Whether this is PASEF data
      @param[in] load_into_memory Whether to load data into memory
      @param[in] irt_trafo_out Output transformation file
      @param[in] irt_mzml_out Output iRT chromatograms file
      
      @return Combined calibration results
    */
    CalibrationResult performLinearThenNonlinearCalibration_(
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
      const String& irt_mzml_out);
      
    /**
      @brief Apply estimated extraction windows to parameters.
      
      Updates the provided extraction parameters with auto-estimated windows
      based on calibration results and user preferences.
      
      @param[in] result Calibration results containing estimated windows
      @param[in,out] cp MS2 extraction parameters to update
      @param[in,out] cp_ms1 MS1 extraction parameters to update
      @param[in] pasef Whether this is PASEF data (for IM applicability)
      @param[in] use_ms1_im Whether MS1 uses ion mobility
    */
    void applyEstimatedWindows_(
      const CalibrationResult& result,
      ChromExtractParams& cp,
      ChromExtractParams& cp_ms1, 
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