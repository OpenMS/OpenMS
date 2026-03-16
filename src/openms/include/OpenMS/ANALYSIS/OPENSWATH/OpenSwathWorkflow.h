// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// Interfaces

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

// Kernel and implementations
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
// #include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
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
      mrm_(false),
      threads_outer_loop_(-1)
    {
    }

    /** @brief Constructor
     *
     *  @param[in] use_ms1_traces Use MS1 data?
     *  @param[in] use_ms1_ion_mobility Use ion mobility extraction on MS1 traces?
     *  @param[out] prm Is data acquired in targeted DIA (e.g. PRM mode) with potentially overlapping windows?
     *  @param[in] pasef Is this diaPASEF data?
     *  @param[in] mrm Is this SRM/MRM data?
     *  @param[in] threads_outer_loop How many threads should be used for the outer
     *  loop (-1 will use all threads in the outer loop)
     *
     *  @note The total number of threads should be divisible by this number
     *  (e.g. use 8 in outer loop if you have 24 threads in total and 3 will be
     *  used for the inner loop).
     *
     *
     **/
    OpenSwathWorkflowBase(bool use_ms1_traces, bool use_ms1_ion_mobility, bool prm, bool pasef, bool mrm, int threads_outer_loop) :
      use_ms1_traces_(use_ms1_traces),
      use_ms1_ion_mobility_(use_ms1_ion_mobility),
      prm_(prm),
      pasef_(pasef),
      mrm_(mrm),
      threads_outer_loop_(threads_outer_loop)
    {
    }

    /** @brief Perform MS1 extraction and store result in ms1_chromatograms
     *
     * @param[in] ms1_map Spectrum Access to the MS1 map
     * @param[in] swath_maps The raw data (swath maps)
     * @param[in] ms1_chromatograms Output vector for MS1 chromatograms
     * @param[in] cp Parameter set for the chromatogram extraction
     * @param[in] transition_exp The set of assays to be extracted and scored
     * @param[in] trafo_inverse Inverse transformation function
     * @param[in] ms1_only If true, will only score on MS1 level and ignore MS2 level
     * @param[in] ms1_isotopes Number of MS1 isotopes to extract (zero means only monoisotopic peak)
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
     * @param[out] chrom_list Output of chromatograms (will be filled with empty chromatogram ptrs)
     * @param[out] coordinates Output of extraction coordinates (will be filled with matching extraction coordinates)
     * @param[out] transition_exp_used The transition experiment used to create the coordinates
     * @param[in] trafo_inverse Inverse transformation function
     * @param[in] cp Parameter set for the chromatogram extraction
     * @param[in] ms1 Whether to perform MS1 (precursor ion) or MS2 (fragment ion) extraction
     * @param[in] ms1_isotopes Number of MS1 isotopes to extract (zero means only monoisotopic peak)     
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
     * case for SRM/MRM / PRM data where often multiple windows with similar (or
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

    /** @brief Whether data is chromatogram-only SRM/MRM data
     *
     * If set to true, indicates that all swath_maps contain only chromatograms
     * (no spectra) and the workflow should use SRM/MRM-specific processing.
    */
    bool mrm_;

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
     *  @param[in] use_ms1_traces Whether to use MS1 data
     *  @param[in] use_ms1_ion_mobility Whether to use ion mobility extraction on MS1 traces
     *  @param[out] prm Whether data is acquired in targeted DIA (e.g. PRM mode) with potentially overlapping windows
     *  @param[in] pasef Is this diaPASEF data?
     *  @param[in] mrm Is this SRM/MRM data?
     *  @param[in] threads_outer_loop How many threads should be used for the outer
     *  loop (-1 will use all threads in the outer loop)
     *
     *  @note The total number of threads should be divisible by this number
     *  (e.g. use 8 in outer loop if you have 24 threads in total and 3 will be
     *  used for the inner loop).
     *
     *
     **/
    OpenSwathWorkflow(bool use_ms1_traces, bool use_ms1_ion_mobility, bool prm, bool pasef, bool mrm, int threads_outer_loop) :
    OpenSwathWorkflowBase(use_ms1_traces, use_ms1_ion_mobility, prm, pasef, mrm, threads_outer_loop)
    {
    }

    /** @brief Load MS1 SpectrumAccessPtr from given swath maps.
     *
     * Searches through the provided swath maps and returns a SpectrumAccessPtr
     * to the first MS1 map found. If no MS1 map is present, returns nullptr.
     *
     * @param[in] swath_maps Vector of SWATH maps to search for MS1 data
     * @param[in] load_into_memory Whether to cache the MS1 map in memory for faster access
     * @return SpectrumAccessPtr to the first MS1 map, or nullptr if no MS1 map exists
     *
     * @note The returned pointer may be cached on disk or in memory depending on load_into_memory parameter
    */
    OpenSwath::SpectrumAccessPtr loadMS1Map(const std::vector<OpenSwath::SwathMap>& swath_maps, bool load_into_memory);

    /** @brief Execute OpenSWATH analysis on a set of SwathMaps and transitions.
     *
     * See OpenSwathWorkflow class for a detailed description of this function.
     *
     * @param[in] swath_maps The raw data (swath maps)
     * @param[in] rt_trafo Retention time transformation description (translating this runs' RT to normalized RT space)
     * @param[in] chromatogram_extraction_params Parameter set for the chromatogram extraction
     * @param[in] ms1_chromatogram_extraction_params Parameter set for the chromatogram extraction of the MS1 data
     * @param[in] feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param[in] assay_library The set of assays to be extracted and scored
     * @param[out] result_featureFile Output feature map to store identified features
     * @param[out] store_features_in_featureFile Whether features should be appended to the output feature map (if this is false, then out_featureFile will be empty)
     * @param[out] result_osw OSW Writer object to store identified features in SQLite format (set store_features to false if using this option)
     * @param[out] result_chromatograms Chromatogram consumer object to store the extracted chromatograms
     * @param[in] batchSize Size of the batches which should be extracted and scored
     * @param[in] ms1_isotopes Number of MS1 isotopes to extract (zero means only monoisotopic peak)
     * @param[in] load_into_memory Whether to cache the current SWATH map in memory
     * @param[in] mrm_mapping_param Parameter for mapping chromatograms to transitions (MRMMapping)
     * @param[in] mobilogram_consumer Optional consumer to write out extracted ion mobilograms
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
               OpenSwathOSWWriter & result_osw,
               Interfaces::IMSDataConsumer * result_chromatograms,
               int batchSize,
               int ms1_isotopes,
               bool load_into_memory,
              const Param & mrm_mapping_param = Param(),
              class MobilogramParquetConsumer * mobilogram_consumer = nullptr);

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
     * @param[out] chromConsumer Chromatogram consumer object to store the extracted chromatograms
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
     *    - Add the identified peak groups to the SQL-based output format (osw_writer)
     *
     * @param[in] ms2_chromatograms Input chromatograms (MS2 level)
     * @param[out] ms1_chromatograms Input chromatograms (MS1-level)
     * @param[in] swath_maps Set of swath map(s) for the current swath window 
     * @param[in] transition_exp The transition experiment (assay library)
     * @param[in] feature_finder_param Parameters for the MRMFeatureFinderScoring
     * @param[in] trafo RT Transformation function
     * @param[in] rt_extraction_window RT extraction window
     * @param[out] output Output map
     * @param[out] osw_writer OSW Writer object to store identified features in SQLite format
     * @param[in] nr_ms1_isotopes Consider this many MS1 isotopes for precursor chromatograms
     * @param[in] ms1only If true, will only score on MS1 level and ignore MS2 level
     * @param[in] mobilogram_consumer Optional consumer to write out extracted ion mobilograms
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
        OpenSwathOSWWriter& osw_writer,
        int nr_ms1_isotopes = 0,
        bool ms1only = false,
        class MobilogramParquetConsumer * mobilogram_consumer = nullptr) const;

    /** @brief Select which compounds to analyze in the next batch (and copy to output)
     *
     * This function will select which compounds or peptides should be analyzed
     * in the current batch (with index "batch_idx"). The selected compounds
     * will be copied into the output structure. The output will contain
     * "batch_size" compounds or peptides.
     *
     * @param[in] transition_exp_used_all The full set of transitions (this will be used to select transitions from)
     * @param[in] transition_exp_used The selected set of transitions (will contain only transitions for the next batch)
     * @param[in] batch_size How many compounds or peptides should be used per batch
     * @param[in] batch_idx Current batch index (only compounds or peptides from batch_idx*batch_size to batch_idx*batch_size+batch_size will be copied)
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
     * @param[in] used_compounds Which peptides or metabolites to be used
     * @param[in] all_transitions Transitions vector from which to select transitions
     * @param[out] output Output vector containing matching transitions (taken from all_transitions)
     *
    */
    void copyBatchTransitions_(const std::vector<OpenSwath::LightCompound>& used_compounds,
      const std::vector<OpenSwath::LightTransition>& all_transitions,
      std::vector<OpenSwath::LightTransition>& output);
  };
}


