// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// Interfaces
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/FORMAT/FileHandler.h> // debug file store only

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
// #include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathTSVWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>

#include <cassert>
#include <limits>

// #define OPENSWATH_WORKFLOW_DEBUG

// The workflow class
namespace OpenMS
{

  /** @brief ChromatogramExtractor parameters
   *
   * A small helper struct to pass the parameters for the chromatogram
   * extraction through to the actual algorithm.
   *
  */
  struct ChromExtractParams
  {
    /// Whether to not extract anything closer than this (in Da) from the upper edge
    double min_upper_edge_dist;
    /// Extraction window in Da or ppm (e.g. 50ppm means extraction +/- 25ppm)
    double mz_extraction_window;
    /// Extraction window in ion mobility
    double im_extraction_window;
    /// Whether the extraction window is given in ppm or Da
    bool ppm;
    /// The extraction function in mass space
    String extraction_function;
    /// The retention time extraction window
    double rt_extraction_window;
    /// Whether to extract some extra in the retention time (can be useful if one wants to look at the chromatogram outside the window)
    double extra_rt_extract;
  };

  class OPENMS_DLLAPI OpenSwathWorkflowBase :
    public ProgressLogger
  {

protected:

    /** @brief Default constructor
     *
     *  Will not use any ms1 traces and use all threads in the outer loop.
     *
     **/
    OpenSwathWorkflowBase() :
      use_ms1_traces_(false),
      use_ms1_ion_mobility_(false),
      prm_(false),
      pasef_(false),
      threads_outer_loop_(-1)
    {
    }

    /** @brief Constructor
     *
     *  @param use_ms1_traces Use MS1 data?
     *  @param use_ms1_ion_mobility Use ion mobility extraction on MS1 traces?
     *  @param prm Is data acquired in targeted DIA (e.g. PRM mode) with potentially overlapping windows?
     *  @param pasef Is this diaPASEF data?
     *  @param threads_outer_loop How many threads should be used for the outer
     *  loop (-1 will use all threads in the outer loop)
     *
     *  @note The total number of threads should be divisible by this number
     *  (e.g. use 8 in outer loop if you have 24 threads in total and 3 will be
     *  used for the inner loop).
     *
     *
     **/
    OpenSwathWorkflowBase(bool use_ms1_traces, bool use_ms1_ion_mobility, bool prm, bool pasef, int threads_outer_loop) :
      use_ms1_traces_(use_ms1_traces),
      use_ms1_ion_mobility_(use_ms1_ion_mobility),
      prm_(prm),
      pasef_(pasef),
      threads_outer_loop_(threads_outer_loop)
    {
    }

    /** @brief Perform MS1 extraction and store result in ms1_chromatograms
     *
     * @param ms1_map Spectrum Access to the MS1 map
     * @param swath_maps The raw data (swath maps)
     * @param ms1_chromatograms Output vector for MS1 chromatograms
     * @param cp Parameter set for the chromatogram extraction
     * @param transition_exp The set of assays to be extracted and scored
     * @param trafo_inverse Inverse transformation function
     * @param ms1_only If true, will only score on MS1 level and ignore MS2 level
     * @param ms1_isotopes Number of MS1 isotopes to extract (zero means only monoisotopic peak)
     *
    */
    void MS1Extraction_(const OpenSwath::SpectrumAccessPtr& ms1_map,
                        const std::vector<OpenSwath::SwathMap>& swath_maps,
                        std::vector<MSChromatogram>& ms1_chromatograms,
                        const ChromExtractParams& cp,
                        const OpenSwath::LightTargetedExperiment& transition_exp,
                        const TransformationDescription& trafo_inverse,
                        bool ms1_only = false,
                        int ms1_isotopes = 0);

    /** @brief Function to prepare extraction coordinates that also correctly handles RT transformations
     *
     * Creates a set of (empty) chromatograms and extraction coordinates with
     * the correct ids, m/z and retention time start/end points to be extracted
     * by the ChromatogramExtractor.
     *
     * Handles RT extraction windows by calculating the correct transformation
     * for each coordinate.
     *
     * @param chrom_list Output of chromatograms (will be filled with empty chromatogram ptrs)
     * @param coordinates Output of extraction coordinates (will be filled with matching extraction coordinates)
     * @param transition_exp_used The transition experiment used to create the coordinates
     * @param trafo_inverse Inverse transformation function
     * @param cp Parameter set for the chromatogram extraction
     * @param ms1 Whether to perform MS1 (precursor ion) or MS2 (fragment ion) extraction
     * @param ms1_isotopes Number of MS1 isotopes to extract (zero means only monoisotopic peak)     
     *
    */
    void prepareExtractionCoordinates_(std::vector< OpenSwath::ChromatogramPtr > & chrom_list,
                                       std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > & coordinates,
                                       const OpenSwath::LightTargetedExperiment & transition_exp_used,
                                       const TransformationDescription& trafo_inverse,
                                       const ChromExtractParams & cp,
                                       const bool ms1 = false,
                                       const int ms1_isotopes = -1) const;


    /**
     * @brief Spectrum Access to the MS1 map (note that this is *not* threadsafe!)
     *
     * @note This pointer is not threadsafe, please use the lightClone() function to create a copy for each thread
     * @note This pointer may be NULL if use_ms1_traces_ is set to false
     *
     */
    OpenSwath::SpectrumAccessPtr ms1_map_ = nullptr;

    /// Whether to use the MS1 traces
    bool use_ms1_traces_;

    /// Whether to use ion mobility extraction on MS1 traces
    bool use_ms1_ion_mobility_;

    /** @brief Whether data is acquired in targeted DIA (e.g. PRM mode) with potentially overlapping windows
     *
     * If set to true, a precursor will only be extracted from a single window
     * that matches in m/z and whose m/z center is *closest* to the library m/z
     * of the precursor. This is required if windows overlap in m/z as is the
     * case for SRM / PRM data where often multiple windows with similar (or
     * overlaping) m/z are used to target different precursors at different RT.
    */
    bool prm_;

    /** @brief Whether data is diaPASEF data
     *
     * If set to true, a precursor will only be extracted from a single window
     * that matches both in m/z and whose ion mobility (drift time) center is
     * *closest* to the library ion mobility of the precursor. This is required
     * if windows overlap in m/z or in ion mobility, as is the case for
     * diaPASEF data.
    */
    bool pasef_;

    /** @brief How many threads should be used for the outer loop
     *
     *  @note A value of -1 will use all threads in the outer loop
     *
     *  @note The total number of threads should be divisible by this number
     *  (e.g. use 8 in outer loop if you have 24 threads in total and 3 will be
     *  used for the inner loop).
     *
     **/
    int threads_outer_loop_;

};

  /**
   * @brief Execute all steps for retention time and m/z calibration of SWATH-MS data
   *
   * Uses a set of robust calibrant peptides (e.g. iRT peptides, common
   * calibrants) perform RT and m/z correction in SWATH-MS data. Currently
   * supports (non-)linear correction of RT against library RT as well
   * as (non-)linear correction of m/z error as a function of m/z.
   *
   * @note The relevant algorithms are implemented in MRMRTNormalizer for RT
   * calibration and SwathMapMassCorrection for m/z calibration.
   *
   * The overall execution flow in this class is as follows (see performRTNormalization() function):
   *   - Extract chromatograms across the whole RT range using simpleExtractChromatograms_()
   *   - Compute calibration functions for RT and m/z using doDataNormalization_()
   *
  */
  class OPENMS_DLLAPI OpenSwathCalibrationWorkflow :
    public OpenSwathWorkflowBase
  {
  public:

    OpenSwathCalibrationWorkflow() :
      OpenSwathWorkflowBase()
    {
    }

    explicit OpenSwathCalibrationWorkflow(bool use_ms1_traces) :
      OpenSwathWorkflowBase(use_ms1_traces, false, false, false, -1)
    {
    }

    /** @brief Perform RT and m/z correction of the input data using RT-normalization peptides.
     *
     * This function extracts the RT normalization chromatograms using
     * simpleExtractChromatograms_() and then uses the chromatograms to find
     * features (in doDataNormalization_()).  If desired, also m/z correction
     * is performed using the lock masses of the given peptides. The provided
     * raw data (swath_maps) are therefore not constant but may be changed in
     * this function.
     *
     * @param irt_transitions A set of transitions used for the RT normalization peptides
     * @param swath_maps The raw data (swath maps)
     * @param[out] im_trafo Ion mobility trafo values on the RT-normalization peptides
     * @param min_rsq Minimal R^2 value that is expected for the RT regression
     * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param cp_irt Parameter set for the chromatogram extraction
     * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param calibration_param Parameter for the m/z and im calibration (see SwathMapMassCorrection)
     * @param debug_level Debug level (writes out the RT normalization chromatograms if larger than 1)
     * @param irt_mzml_out Output Chromatogram mzML containing the iRT peptides (if not empty,
     *        iRT chromatograms will be stored in this file)
     * @param pasef whether the data is PASEF data (should match transitions by their IM)
     * @param load_into_memory Whether to cache the current SWATH map in memory
     *
    */
    TransformationDescription performRTNormalization(const OpenSwath::LightTargetedExperiment & irt_transitions,
      std::vector< OpenSwath::SwathMap > & swath_maps,
      TransformationDescription& im_trafo,
      double min_rsq,
      double min_coverage,
      const Param & feature_finder_param,
      const ChromExtractParams & cp_irt,
      const Param& irt_detection_param,
      const Param& calibration_param,
      const String& irt_mzml_out,
      Size debug_level,
      bool pasef = false,
      bool load_into_memory = false);

  public:

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
     * @param transition_exp_ The transitions for the normalization peptides
     * @param chromatograms The extracted chromatograms
     * @param[out] im_trafo Ion mobility trafo values on the RT-normalization peptides
     * @param swath_maps The raw data (swath maps)     
     * @param min_rsq Minimal R^2 value that is expected for the RT regression
     * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param default_ffparam Parameter set for the feature finding in chromatographic dimension
     * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param calibration_param Parameter for the m/z and im calibration (see SwathMapMassCorrection)
     * @param pasef whether this data is pasef data with potentially overlapping m/z windows (differing by IM)
     *
     * @note This function is based on the algorithm inside the OpenSwathRTNormalizer tool
     *
    */
    TransformationDescription doDataNormalization_(const OpenSwath::LightTargetedExperiment& transition_exp_,
      const std::vector< OpenMS::MSChromatogram >& chromatograms,
      TransformationDescription& im_trafo,
      std::vector< OpenSwath::SwathMap > & swath_maps,
      double min_rsq,
      double min_coverage,
      const Param& default_ffparam,
      const Param& irt_detection_param,
      const Param& calibration_param,
      const bool pasef);

    /** @brief Simple method to extract chromatograms (for the RT-normalization peptides)
     *
     * @param swath_maps The raw data (swath maps)
     * @param irt_transitions A set of transitions used for the RT normalization peptides
     * @param chromatograms The extracted chromatograms (output)
     * @param trafo Transformation description for RT normalization
     * @param cp Parameter set for the chromatogram extraction
     * @param load_into_memory Whether to cache the current SWATH map in memory
     * @param pasef whether the data is PASEF data with possible overlapping m/z windows (with different ion mobility)
     *
    */
    void simpleExtractChromatograms_(const std::vector< OpenSwath::SwathMap > & swath_maps,
                                     const OpenSwath::LightTargetedExperiment & irt_transitions,
                                     std::vector< OpenMS::MSChromatogram > & chromatograms,
                                     const TransformationDescription& trafo,
                                     const ChromExtractParams & cp,
                                     bool pasef,
                                     bool load_into_memory);

    /** @brief Add two chromatograms
     *
     * @param base_chrom The base chromatogram to which we will add intensity
     * @param newchrom The chromatogram to be added
     *
    */
    static void addChromatograms(MSChromatogram& base_chrom, const MSChromatogram& newchrom);

  };

  /**
   * @brief Execute all steps in an \ref TOPP_OpenSwathWorkflow "OpenSwath" analysis
   *
   * The workflow will perform a complete OpenSWATH analysis. Optionally,
   * a calibration of m/z and retention time (mapping peptides to normalized
   * space and correcting m/z error) can be performed beforehand using the
   * OpenSwathCalibrationWorkflow class.
   *
   * For diaPASEF workflows where ion mobility windows are overlapping, precursors may be found in multiple SWATHs.
   * In this case, precursors are only extracted from the SWATH in which they are most centered across ion mobility
   * (Provided -pasef flag is set).
   *
   * The overall execution flow in this class is as follows (see performExtraction() function)
   *
   *    - Obtain precursor ion chromatograms (if enabled) through MS1Extraction_()
   *    - Perform scoring of precursor ion chromatograms if no MS2 is given
   *    - Iterate through each SWATH-MS window:
   *      - Select which transitions to extract (proceed in batches) using OpenSwathHelper::selectSwathTransitions()
   *      - Iterate through each batch of transitions:
   *        - Extract current batch of transitions from current SWATH window:
   *          - Select transitions for current batch (see selectCompoundsForBatch_())
   *          - Prepare transition extraction (see prepareExtractionCoordinates_())
   *          - Extract transitions using ChromatogramExtractor::extractChromatograms()
   *          - Convert data to OpenMS format using ChromatogramExtractor::return_chromatogram()
   *        - Score extracted transitions (see scoreAllChromatograms_())
   *        - Write scored chromatograms and peak groups to disk (see writeOutFeaturesAndChroms_())
   *
   */
  class OPENMS_DLLAPI OpenSwathWorkflow :
    public OpenSwathWorkflowBase
  {
    typedef OpenSwath::LightTransition TransitionType;
    typedef MRMTransitionGroup< MSChromatogram, TransitionType> MRMTransitionGroupType;

  public:

    /** @brief Constructor
     *
     *  @param use_ms1_traces Whether to use MS1 data
     *  @param use_ms1_ion_mobility Whether to use ion mobility extraction on MS1 traces
     *  @param prm Whether data is acquired in targeted DIA (e.g. PRM mode) with potentially overlapping windows
     *  @param pasef Is this diaPASEF data?
     *  @param threads_outer_loop How many threads should be used for the outer
     *  loop (-1 will use all threads in the outer loop)
     *
     *  @note The total number of threads should be divisible by this number
     *  (e.g. use 8 in outer loop if you have 24 threads in total and 3 will be
     *  used for the inner loop).
     *
     *
     **/
    OpenSwathWorkflow(bool use_ms1_traces, bool use_ms1_ion_mobility, bool prm, bool pasef, int threads_outer_loop) :
    OpenSwathWorkflowBase(use_ms1_traces, use_ms1_ion_mobility, prm, pasef, threads_outer_loop)
    {
    }

    /** @brief Execute OpenSWATH analysis on a set of SwathMaps and transitions.
     *
     * See OpenSwathWorkflow class for a detailed description of this function.
     *
     * @param swath_maps The raw data (swath maps)
     * @param rt_trafo Retention time transformation description (translating this runs' RT to normalized RT space)
     * @param chromatogram_extraction_params Parameter set for the chromatogram extraction
     * @param ms1_chromatogram_extraction_params Parameter set for the chromatogram extraction of the MS1 data
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param assay_library The set of assays to be extracted and scored
     * @param result_featureFile Output feature map to store identified features
     * @param store_features_in_featureFile Whether features should be appended to the output feature map (if this is false, then out_featureFile will be empty)
     * @param result_tsv TSV Writer object to store identified features in csv format (set store_features to false if using this option)
     * @param result_osw OSW Writer object to store identified features in SQLite format (set store_features to false if using this option)
     * @param result_chromatograms Chromatogram consumer object to store the extracted chromatograms
     * @param batchSize Size of the batches which should be extracted and scored
     * @param ms1_isotopes Number of MS1 isotopes to extract (zero means only monoisotopic peak)
     * @param load_into_memory Whether to cache the current SWATH map in memory
     *
     * @note Speed and memory performance can be influenced by \p batchSize and
     * \p load_into_memory where larger batch sizes increase memory and
     * potentially decrease the utility of parallelization while loading data
     * into memory will increase memory usage but decrease execution time.
     *
    */
    void performExtraction(const std::vector<OpenSwath::SwathMap>& swath_maps,
                           const TransformationDescription& rt_trafo,
                           const ChromExtractParams & chromatogram_extraction_params,
                           const ChromExtractParams & ms1_chromatogram_extraction_params,
                           const Param & feature_finder_param,
                           const OpenSwath::LightTargetedExperiment& assay_library,
                           FeatureMap& result_featureFile,
                           bool store_features_in_featureFile,
                           OpenSwathTSVWriter & result_tsv,
                           OpenSwathOSWWriter & result_osw,
                           Interfaces::IMSDataConsumer * result_chromatograms,
                           int batchSize,
                           int ms1_isotopes,
                           bool load_into_memory);

  protected:


    /** @brief Write output features and chromatograms
     *
     * Writes output chromatograms to the provided chromatogram consumer
     * (presumably to disk) and output features to the provided FeatureMap.
     *
     * @param[in] chromatograms Output chromatograms to be passed to the consumer
     * @param[in] ms1_chromatograms Output chromatograms (MS1 level) to be passed to the consumer
     * @param[in] featureFile Features to be appended to the @p out_featureFile
     * @param[out] out_featureFile Output FeatureMap to which the features will be appended
     * @param[in] store_features Whether features should be appended to the output
     *        feature map (if this is false, then out_featureFile will be empty)
     * @param chromConsumer Chromatogram consumer object to store the extracted chromatograms
     *
     * @note This should be wrapped in an OpenMP critical block
    */
    void writeOutFeaturesAndChroms_(std::vector< OpenMS::MSChromatogram > & chromatograms,
                                    std::vector< MSChromatogram >& ms1_chromatograms,
                                    const FeatureMap & featureFile,
                                    FeatureMap& out_featureFile,
                                    bool store_features,
                                    Interfaces::IMSDataConsumer * chromConsumer);

    /** @brief Perform scoring on a set of chromatograms
     *
     *  This will generate a new object of type MRMTransitionGroup for each
     *  compound or peptide in the provided assay library and link the
     *  transition meta information with the extracted chromatograms. This will
     *  then be used to perform peak picking and peak scoring through
     *  MRMTransitionGroupPicker and MRMFeatureFinderScoring. The assay library
     *  is provided as transition_exp and the chromatograms in
     *  ms2_chromatograms.
     *
     * The overall execution flow is as follows:
     *
     *  - Iterate over all assays (compounds / peptides) in the transition_exp
     *    - Create a new MRMTransitionGroup
     *    - Iterate over all transitions in an assay
     *      - Find the relevant chromatogram for the given transition, convert it and filter it by RT
     *      - Add the chromatogram and transition to the MRMTransitionGroup
     *    - Add a single MS1 chromatogram of the mono-isotopic precursor to the
     *    MRMTransitionGroup, if available (named "groupId_Precursor_i0")
     *    - Find peakgroups in the chromatogram set (see MRMTransitionGroupPicker::pickTransitionGroup)
     *    - Score peakgroups in the chromatogram set (see MRMFeatureFinderScoring::scorePeakgroups)
     *    - Add the identified peak groups to the TSV writer (tsv_writer) and the SQL-based output format (osw_writer)
     *
     * @param ms2_chromatograms Input chromatograms (MS2 level)
     * @param ms1_chromatograms Input chromatograms (MS1-level)
     * @param swath_maps Set of swath map(s) for the current swath window 
     * @param transition_exp The transition experiment (assay library)
     * @param feature_finder_param Parameters for the MRMFeatureFinderScoring
     * @param trafo RT Transformation function
     * @param rt_extraction_window RT extraction window
     * @param output Output map
     * @param tsv_writer TSV writer for storing output (on the fly)
     * @param osw_writer OSW Writer object to store identified features in SQLite format
     * @param nr_ms1_isotopes Consider this many MS1 isotopes for precursor chromatograms
     * @param ms1only If true, will only score on MS1 level and ignore MS2 level
     *
    */
    void scoreAllChromatograms_(
        const std::vector<OpenMS::MSChromatogram>& ms2_chromatograms,
        const std::vector<OpenMS::MSChromatogram>& ms1_chromatograms,
        const std::vector<OpenSwath::SwathMap>& swath_maps,
        const OpenSwath::LightTargetedExperiment& transition_exp,
        const Param& feature_finder_param,
        const TransformationDescription& trafo,
        const double rt_extraction_window,
        FeatureMap& output,
        OpenSwathTSVWriter& tsv_writer,
        OpenSwathOSWWriter& osw_writer,
        int nr_ms1_isotopes = 0,
        bool ms1only = false) const;

    /** @brief Select which compounds to analyze in the next batch (and copy to output)
     *
     * This function will select which compounds or peptides should be analyzed
     * in the current batch (with index "batch_idx"). The selected compounds
     * will be copied into the output structure. The output will contain
     * "batch_size" compounds or peptides.
     *
     * @param transition_exp_used_all The full set of transitions (this will be used to select transitions from)
     * @param transition_exp_used The selected set of transitions (will contain only transitions for the next batch)
     * @param batch_size How many compounds or peptides should be used per batch
     * @param batch_idx Current batch index (only compounds or peptides from batch_idx*batch_size to batch_idx*batch_size+batch_size will be copied)
     *
     * @note The proteins will be copied completely without checking for a match
     *
    */
    void selectCompoundsForBatch_(const OpenSwath::LightTargetedExperiment& transition_exp_used_all,
      OpenSwath::LightTargetedExperiment& transition_exp_used, int batch_size, size_t batch_idx);

    /** @brief Helper function for selectCompoundsForBatch_()
     *
     * Copy all transitions matching to one of the compounds in the selected
     * peptide vector from all_transitions to the output.
     *
     * @param used_compounds Which peptides or metabolites to be used
     * @param all_transitions Transitions vector from which to select transitions
     * @param output Output vector containing matching transitions (taken from all_transitions)
     *
    */
    void copyBatchTransitions_(const std::vector<OpenSwath::LightCompound>& used_compounds,
      const std::vector<OpenSwath::LightTransition>& all_transitions,
      std::vector<OpenSwath::LightTransition>& output);
  };
}


