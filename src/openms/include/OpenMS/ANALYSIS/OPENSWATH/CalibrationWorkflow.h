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
    
    - NULL_TRANSFORMATION: No iRT data available or explicit skip (represents a null/no-transformation case)
    - STATIC_FILES: Use pre-defined iRT peptide libraries from files (same for all runs). If all samples have spiked in iRTs, or common expected peptides.
    - SAMPLE_ONCE: Sample iRT peptides once from transition library, reuse for all runs  
    - SAMPLE_PER_RUN: Sample fresh iRT peptides independently for each run
    - RUN_SPECIFIC: Use different iRT files per run (each run has its own iRT files). These are usually generated from accompanying pseudo-DDA spectra analysis
  */
  enum class IrtStrategy 
  {
    NULL_TRANSFORMATION,
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
  */
  class OPENMS_DLLAPI CalibrationWorkflow : 
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    /// Results from calibration workflow
    struct CalibrationResult 
    {
      // Transformation
      /// RT normalization transformation (fitted)
      TransformationDescription rt_trafo;
      /// Ion mobility transformation (fitted). May be empty if no IM calibration performed.
      TransformationDescription im_trafo;
      
      // Estimated Extraction Windows
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
      2. If full transition library is provided and auto-iRT is enabled → SAMPLE_ONCE or SAMPLE_PER_RUN
      3. Otherwise → NULL_TRANSFORMATION (no iRT data available)
      
      @param[in] full_transition_exp Full transition experiment for auto-sampling (empty if not available)
      @param[in] num_runs Total number of runs to process (affects strategy choice)
      @return The determined IRT strategy
      
      @note When NULL_TRANSFORMATION is returned, callers can decide whether to skip calibration
      or throw an exception based on their requirements
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
      @param[in] debug_level Debug level (0 = no debug output, >1 = verbose)
      
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
      const String& irt_mzml_out = "",
      Size debug_level = 0
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

    /** @brief Perform RT and m/z correction of the input data using RT-normalization peptides.
     *
     * This function extracts the RT normalization chromatograms using
     * simpleExtractChromatograms_() and then uses the chromatograms to find
     * features (in doDataNormalization_()).  If desired, also m/z correction
     * is performed using the lock masses of the given peptides. The provided
     * raw data (swath_maps) are therefore not constant but may be changed in
     * this function.
     *
     * @param[in] irt_transitions A set of transitions used for the RT normalization peptides
     * @param[in] swath_maps The raw data (swath maps)
     * @param[out] im_trafo Ion mobility trafo values on the RT-normalization peptides
     * @param[in] min_rsq Minimal R^2 value that is expected for the RT regression
     * @param[in] min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param[in] feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param[in] cp_irt Parameter set for the chromatogram extraction
     * @param[in] irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param[in] calibration_param Parameter for the m/z and im calibration (see SwathMapMassCorrection)
     * @param[in] mrm_mapping_param Parameter for mapping chromatograms to transitions (MRMMapping)
     * @param[in] debug_level Debug level (writes out the RT normalization chromatograms if larger than 1)
     * @param[out] irt_mzml_out Output Chromatogram mzML containing the iRT peptides (if not empty,
     *        iRT chromatograms will be stored in this file)
     * @param[in] pasef whether the data is PASEF data (should match transitions by their IM)
     * @param[in] load_into_memory Whether to cache the current SWATH map in memory
     *
    */
    TransformationDescription performRTNormalization(
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
      bool pasef = false,
      bool load_into_memory = false);

    /** @brief Perform retention time and m/z calibration
     *
     * Uses MRMRTNormalizer for RT calibration and SwathMapMassCorrection for m/z calibration.
     *
     * The overall execution flow is as follows:
     *   - Estimate the retention time range of the iRT peptides over all assays (see OpenSwathHelper::estimateRTRange())
     *   - Store the peptide retention times in an intermediate map
     *   - Pick input chromatograms to identify RT pairs from the input data
     *   using MRMFeatureFinderScoring, which will be used without the RT
     *   scoring enabled
     *   - Find most likely correct feature for each compound (see OpenSwathHelper::simpleFindBestFeature())
     *   - Perform the outlier detection (see MRMRTNormalizer)
     *   - Check whether the found peptides fulfill the binned coverage criteria set by the user.
     *   - Select the "correct" peaks for m/z correction (e.g. remove those not
     *   part of the linear regression)
     *   - Perform m/z and IM calibration (see SwathMapMassCorrection)
     *   - Store transformation, using the selected model
     *
     * @param[in] targeted_exp The transitions for the normalization peptides
     * @param[out] chromatograms The extracted chromatograms
     * @param[out] im_trafo Ion mobility trafo values on the RT-normalization peptides
     * @param[in] swath_maps The raw data (swath maps)     
     * @param[in] min_rsq Minimal R^2 value that is expected for the RT regression
     * @param[in] min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param[in] default_ffparam Parameter set for the feature finding in chromatographic dimension
     * @param[in] irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param[in] calibration_param Parameter for the m/z and im calibration (see SwathMapMassCorrection)
     * @param[in] pasef whether this data is pasef data with potentially overlapping m/z windows (differing by IM)
     *
     * @note This function is based on the algorithm inside the OpenSwathRTNormalizer tool
     *
    */
    TransformationDescription doDataNormalization_(
      const OpenSwath::LightTargetedExperiment& targeted_exp,
      const std::vector< OpenMS::MSChromatogram >& chromatograms,
      TransformationDescription& im_trafo,
      std::vector< OpenSwath::SwathMap > & swath_maps,
      double min_rsq,
      double min_coverage,
      const Param& default_ffparam,
      const Param& irt_detection_param,
      const Param& calibration_param,
      const bool pasef);

    /** @brief Get estimated MS2 m/z extraction window.
     *
     * Returns the MS2 m/z extraction window (full width, in ppm) estimated
     * during calibration. This value is only valid after performCalibration()
     * has been called and window estimation was enabled.
     *
     * @return Estimated MS2 m/z window in ppm, or -1.0 if not computed
    */
    double getEstimatedMzWindow() const { return estimated_mz_window_; }

    /** @brief Get estimated MS2 ion mobility extraction window.
     *
     * Returns the MS2 ion mobility extraction window (full width) estimated
     * during calibration. This value is only valid after performCalibration()
     * has been called on PASEF/ion mobility data with window estimation enabled.
     *
     * @return Estimated MS2 IM window, or -1.0 if not applicable/computed
    */
    double getEstimatedImWindow() const { return estimated_im_window_; }

    /** @brief Get estimated MS1 m/z extraction window.
     *
     * Returns the MS1 m/z extraction window (full width, in ppm) estimated
     * during calibration. This value is only valid after performCalibration()
     * has been called and MS1 window estimation was enabled.
     *
     * @return Estimated MS1 m/z window in ppm, or -1.0 if not computed
    */
    double getEstimatedMs1MzWindow() const { return estimated_ms1_mz_window_; }

    /** @brief Get estimated MS1 ion mobility extraction window.
     *
     * Returns the MS1 ion mobility extraction window (full width) estimated
     * during calibration. This value is only valid after performCalibration()
     * has been called on PASEF/ion mobility data with MS1 IM window estimation enabled.
     *
     * @return Estimated MS1 IM window, or -1.0 if not applicable/computed
    */
    double getEstimatedMs1ImWindow() const { return estimated_ms1_im_window_; }

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

    // Estimated windows stored during calibration (internal caching)
    double estimated_mz_window_{-1.0};
    double estimated_im_window_{-1.0};
    double estimated_ms1_mz_window_{-1.0};
    double estimated_ms1_im_window_{-1.0};
    
    /// @name Private member variables for parameter caching  
    //@{
    
    // iRT File Parameters
    String linear_irt_file_;
    String nonlinear_irt_file_;
    // Run-specific iRT file lists (positional mapping: nth entry -> nth run)
    std::vector<String> linear_irt_files_list_;
    std::vector<String> nonlinear_irt_files_list_;
    
    // Auto-iRT Sampling Parameters
    bool auto_irt_enabled_;
    int auto_irt_irt_bins_;
    int auto_irt_irt_peptides_per_bin_;
    int auto_irt_irt_seed_;
    int auto_irt_irt_bins_nonlinear_;
    int auto_irt_irt_peptides_per_bin_nonlinear_;
    double auto_irt_linear_top_fraction_;
    double auto_irt_nonlinear_top_fraction_;
    
    // Linear Calibration Parameters
    String linear_outlier_detection_;
      
    // Nonlinear Calibration Parameters
    String nonlinear_outlier_detection_;
    
    // Window Estimation Parameters
    bool windows_estimate_rt_;
    bool windows_estimate_mz_;
    bool windows_estimate_im_;
    double windows_rt_percentile_;
    double rt_estimation_padding_factor_;
    
    // Quality Control Parameters
    double min_rsq_;
    double min_coverage_;
    
    //@}
    
    /// @name Private implementation methods
    //@{
    
    /**
      @brief Perform linear-only calibration workflow.
      
      @param[in,out] swath_maps SWATH data maps
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
      @param[in] debug_level Debug level (0 = no debug output, >1 = verbose)
      
      @return Linear calibration results
    */
    CalibrationResult performLinearCalibration_(
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
      Size debug_level);
    
    /**
      @brief Perform linear + nonlinear calibration workflow.
      
      First performs a linear calibration, then applies nonlinear refinement
      using a separate set of iRT transitions.
      
      @param[in,out] swath_maps SWATH data maps
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
      @param[in] debug_level Debug level (0 = no debug output, >1 = verbose)
      
      @return Combined calibration results
    */
    CalibrationResult performLinearThenNonlinearCalibration_(
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
      Size debug_level);
      
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

      Behavior:
        - If @p applicable is false (e.g., no IM data), logs an INFO and leaves @p dst_param unchanged.
        - If the estimate is invalid, logs a WARN and leaves @p dst_param unchanged.
        - If the estimate is valid and @p commit is true, logs an INFO and assigns @p dst_param = @p estimate.
        - If the estimate is valid and @p commit is false, logs an INFO that reports the estimate and that the user value is kept.

      Typical usage:
        - RT window (seconds)
        - MS2 m/z window (ppm)
        - MS1 m/z window (ppm)
        - IM window (1/k0), only when applicable (e.g., PASEF/IM data)
      
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
      @brief Check if an estimated extraction window value is valid.

      A window is considered valid if it is finite and strictly greater than a
      small positive threshold. This guards against denormals (e.g., ~1e-310),
      zeros, negative values, and NaNs/Inf.
      
      @param[in] v Window value to check
      @param[in] min_positive Minimum positive threshold
      @return True if window is valid for use
    */
    inline bool isValidWindow_(double v, double min_positive = 1e-9) const noexcept;
    //@}
  };

} // namespace OpenMS