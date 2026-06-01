// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#define USE_SP_INTERFACE

// Actual scoring
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/FEATUREFINDER/EmgScoring.h>

// Kernel classes
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
  // Forward declaration for optional mobilogram consumer
  class MobilogramParquetConsumer;
  class ProteaseDigestion;


  /**
  @brief The MRMFeatureFinder finds and scores peaks of transitions that co-elute.

  It does so using an internal peakpicker for each chromatogram and then
  creating consensus / meta-peaks (MRMFeatures) that contain the information of
  all corresponding chromatograms at the peak-position. It then goes on to
  score those MRMFeatures using different criteria described in the
  MRMScoring class.

  Internally, all peak group detection is performed in MRMTransitionGroupPicker
  which segments the data and determines consensus peaks across traces
  (MRMFeatures). All scoring is delegated to the OpenSwathScoring class which
  implements i) chromatographic scores, ii) library based scores and iii) full
  spectrum (DIA) scores. These scores are retrieved from the OpenSwathScoring
  class and added to the MRMFeatures found in this algorithm. Note that the
  OpenSwathScoring is a facade that can be used to communicate with the
  underlying actual scoring engines and assembles the scores inside a scoring
  object called OpenSwath_Scores where they are easy to retrieve.

  @htmlinclude OpenMS_MRMFeatureFinderScoring.parameters

  */
  class OPENMS_DLLAPI MRMFeatureFinderScoring :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:
    ///Type definitions
    //@{
    typedef OpenSwath::LightTransition TransitionType;
    typedef OpenSwath::LightTargetedExperiment TargetedExpType;
    typedef OpenSwath::LightCompound PeptideType;
    typedef OpenSwath::LightProtein ProteinType;
    typedef OpenSwath::LightModification ModificationType;
    // a transition group holds the chromatographic data and peaks across
    // multiple chromatograms from the same compound
    typedef MRMTransitionGroup< MSChromatogram, TransitionType> MRMTransitionGroupType;
    typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;

    //@}

    /// Constructor
    MRMFeatureFinderScoring();

    /// Destructor
    ~MRMFeatureFinderScoring() override;

    /// Picker and prepare functions
    //@{

    /** @brief Pick and score features in a single experiment from chromatograms
     *
     * Function for wrapping in Python, only uses OpenMS datastructures and
     * does not return the map.
     *
     * @param[in] chromatograms The input chromatograms
     * @param[out] output The output features with corresponding scores
     * @param[in] transition_exp The transition list describing the experiment
     * @param[in] trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the transition list
     * @param[in] swath_map Optional SWATH-MS (DIA) map corresponding from which the chromatograms were extracted
     *
    */
    void pickExperiment(const PeakMap & chromatograms,
                        FeatureMap& output,
                        const TargetedExperiment& transition_exp,
                        const TransformationDescription& trafo,
                        const PeakMap& swath_map);

    /** @brief Pick and score features in a single experiment from chromatograms
     *
     * @param[in] input The input chromatograms
     * @param[out] output The output features with corresponding scores
     * @param[in] transition_exp The transition list describing the experiment
     * @param[in] trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param[in] swath_maps Optional SWATH-MS (DIA) map corresponding from which
     *                  the chromatograms were extracted. Use empty map if no
     *                  data is available.
     * @param[in] transition_group_map Output mapping of transition groups
     *
    */
    void pickExperiment(const OpenSwath::SpectrumAccessPtr& input,
                        FeatureMap& output,
                        const OpenSwath::LightTargetedExperiment& transition_exp,
                        const TransformationDescription& trafo,
                        const std::vector<OpenSwath::SwathMap>& swath_maps,
                        TransitionGroupMapType& transition_group_map);

    /** @brief Prepares the internal mappings of peptides and proteins.
     *
     * Calling this method _is_ required before calling scorePeakgroups.
     *
     * @param[in] transition_exp The transition list describing the experiment
     *
    */
    void prepareProteinPeptideMaps_(const OpenSwath::LightTargetedExperiment& transition_exp);
    //@}

    /** @brief Score all peak groups of a transition group
     *
     * Iterate through all features found along the chromatograms of the
     * transition group and score each one individually.
     *
     * @param[in] transition_group The MRMTransitionGroup to be scored (input)
     * @param[in] trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param[in] swath_maps Optional SWATH-MS (DIA) map corresponding from which
     *                   the chromatograms were extracted. Use empty map if no
     *                   data is available.
     * @param[out] output The output features with corresponding scores (the found
     *               features will be added to this FeatureMap).
     * @param[in] ms1only Whether to only do MS1 scoring and skip all MS2 scoring
     * @param[in] mobilogram_consumer Optional consumer to write out extracted ion mobilograms
     *
    */
  void scorePeakgroups(MRMTransitionGroupType& transition_group,
             const TransformationDescription & trafo,
             const std::vector<OpenSwath::SwathMap>& swath_maps,
             FeatureMap& output,
             bool ms1only = false,
             MobilogramParquetConsumer* mobilogram_consumer = nullptr) const;

    /** @brief Set the flag for strict mapping
    */
    void setStrictFlag(bool f)
    {
      strict_ = f;
    }

    /** @brief Add an MS1 map containing spectra
     *
     * For DIA (SWATH-MS), an optional MS1 map can be supplied which can be
     * used to extract precursor ion signal and provides additional scores. If
     * no MS1 map is provided, the respective scores are not calculated.
     *
     * @param[in] ms1_map The raw mass spectrometric MS1 data
     *
    */
    void setMS1Map(OpenSwath::SpectrumAccessPtr ms1_map)
    {
      ms1_map_ = ms1_map;
    }

    /** @brief Map the chromatograms to the transitions.
     *
     * Map an input chromatogram experiment (mzML) and transition list (TraML)
     * onto each other when they share identifiers, e.g. if the transition id
     * is the same as the chromatogram native id.
     *
     * @param[in] input The input chromatograms
     * @param[in] transition_exp The transition list describing the experiment
     * @param[in] transition_group_map Mapping of transition groups
     * @param[in] trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param[in] rt_extraction_window The used retention time extraction window
     *
    */
    void mapExperimentToTransitionList(const OpenSwath::SpectrumAccessPtr& input,
                                       const OpenSwath::LightTargetedExperiment& transition_exp,
                                       TransitionGroupMapType& transition_group_map,
                                       TransformationDescription trafo,
                                       double rt_extraction_window);
private:

    /**
     * @brief Cache of transition-group invariant data to avoid per-feature recomputation
     */
    struct TransitionGroupCache
    {
      std::vector<double> normalized_library_intensity;
      std::vector<std::string> transition_native_ids;
      std::vector<std::string> precursor_ids;
      std::vector<String> transition_native_ids_openms;
      std::vector<String> precursor_ids_openms;
    };

    /**
     * @brief Pooled variant of OpenSwath_Ind_Scores to allow pre-allocation and reuse
     */
    struct OpenSwath_Ind_Scores_Pooled : public OpenSwath_Ind_Scores
    {
      /// Reserve storage for all identification score vectors.
      void preallocate(Size capacity)
      {
        ind_transition_names.reserve(capacity);
        ind_isotope_correlation.reserve(capacity);
        ind_isotope_overlap.reserve(capacity);
        ind_massdev_score.reserve(capacity);
        ind_xcorr_coelution_score.reserve(capacity);
        ind_xcorr_shape_score.reserve(capacity);
        ind_log_sn_score.reserve(capacity);
        ind_area_intensity.reserve(capacity);
        ind_total_area_intensity.reserve(capacity);
        ind_intensity_score.reserve(capacity);
        ind_apex_intensity.reserve(capacity);
        ind_apex_position.reserve(capacity);
        ind_fwhm.reserve(capacity);
        ind_total_mi.reserve(capacity);
        ind_log_intensity.reserve(capacity);
        ind_intensity_ratio.reserve(capacity);
        ind_mi_ratio.reserve(capacity);
        ind_mi_score.reserve(capacity);
        ind_im_drift.reserve(capacity);
        ind_im_drift_left.reserve(capacity);
        ind_im_drift_right.reserve(capacity);
        ind_im_delta.reserve(capacity);
        ind_im_delta_score.reserve(capacity);
        ind_im_log_intensity.reserve(capacity);
        ind_im_contrast_coelution.reserve(capacity);
        ind_im_contrast_shape.reserve(capacity);
        ind_im_sum_contrast_coelution.reserve(capacity);
        ind_im_sum_contrast_shape.reserve(capacity);
        ind_start_position_at_5.reserve(capacity);
        ind_end_position_at_5.reserve(capacity);
        ind_start_position_at_10.reserve(capacity);
        ind_end_position_at_10.reserve(capacity);
        ind_start_position_at_50.reserve(capacity);
        ind_end_position_at_50.reserve(capacity);
        ind_total_width.reserve(capacity);
        ind_tailing_factor.reserve(capacity);
        ind_asymmetry_factor.reserve(capacity);
        ind_slope_of_baseline.reserve(capacity);
        ind_baseline_delta_2_height.reserve(capacity);
        ind_points_across_baseline.reserve(capacity);
        ind_points_across_half_height.reserve(capacity);
      }

      /// Clear all identification scores while preserving allocated capacity.
      void reset()
      {
        ind_num_transitions = 0;
        ind_transition_names.clear();
        ind_isotope_correlation.clear();
        ind_isotope_overlap.clear();
        ind_massdev_score.clear();
        ind_xcorr_coelution_score.clear();
        ind_xcorr_shape_score.clear();
        ind_log_sn_score.clear();
        ind_area_intensity.clear();
        ind_total_area_intensity.clear();
        ind_intensity_score.clear();
        ind_apex_intensity.clear();
        ind_apex_position.clear();
        ind_fwhm.clear();
        ind_total_mi.clear();
        ind_log_intensity.clear();
        ind_intensity_ratio.clear();
        ind_mi_ratio.clear();
        ind_mi_score.clear();
        ind_im_drift.clear();
        ind_im_drift_left.clear();
        ind_im_drift_right.clear();
        ind_im_delta.clear();
        ind_im_delta_score.clear();
        ind_im_log_intensity.clear();
        ind_im_contrast_coelution.clear();
        ind_im_contrast_shape.clear();
        ind_im_sum_contrast_coelution.clear();
        ind_im_sum_contrast_shape.clear();
        ind_start_position_at_5.clear();
        ind_end_position_at_5.clear();
        ind_start_position_at_10.clear();
        ind_end_position_at_10.clear();
        ind_start_position_at_50.clear();
        ind_end_position_at_50.clear();
        ind_total_width.clear();
        ind_tailing_factor.clear();
        ind_asymmetry_factor.clear();
        ind_slope_of_baseline.clear();
        ind_baseline_delta_2_height.clear();
        ind_points_across_baseline.clear();
        ind_points_across_half_height.clear();
      }
    };

    /** @brief Splits combined transition groups into detection transition groups
     *
     * For standard assays, transition_group_detection is identical to transition_group and the others are empty.
     *
     * @param[in] transition_group Containing all detecting, identifying transitions
     * @param[out] transition_group_detection To be filled with detecting transitions
    */
    void splitTransitionGroupsDetection_(const MRMTransitionGroupType& transition_group,
                                         MRMTransitionGroupType& transition_group_detection) const;

    /** @brief Splits combined transition groups into identification transition groups
     *
     * For standard assays, transition_group_identification is empty. When UIS scoring
     * is enabled, it contains the corresponding identification transitions.
     *
     * @param[in] transition_group Containing all detecting, identifying transitions
     * @param[out] transition_group_identification To be filled with identifying transitions
     * @param[out] transition_group_identification_decoy To be filled with identifying decoy transitions
    */
    void splitTransitionGroupsIdentification_(const MRMTransitionGroupType& transition_group,
                                              MRMTransitionGroupType& transition_group_identification,
                                              MRMTransitionGroupType& transition_group_identification_decoy) const;

    /** @brief Provides scoring for target and decoy identification against detecting transitions
     *
     * The function is used twice, for target and decoy identification transitions. The results are
     * reported analogously to the ones for detecting transitions but must be stored separately.
     *
     * @param[in,out] transition_group_identification Containing all detecting and identifying transitions
     * @param[in,out] transition_group_detection Containing all detecting transitions
     * @param[in] scorer An instance of OpenSwathScoring
     * @param[in] feature_idx The index of the current feature
     * @param[in] feature_id The id of the current feature
     * @param[in] native_ids_detection The native IDs of the detecting transitions
     * @param[in] signal_noise_estimators_identification Precomputed signal-to-noise estimators for the identification transitions
     * @param[in] det_intensity_ratio_score The intensity score of the detection transitions for normalization
     * @param[in] det_mi_ratio_score The MI score of the detection transitions for normalization
     * @param[in] swath_maps Optional SWATH-MS (DIA) map corresponding from which
     *                  the chromatograms were extracted. Use empty map if no
     *                  data is available.
     * @param[in] drift_target The target drift value
     * @param[in] im_range Ion mobility subrange to consider (used as filter); can be empty (i.e. no IM filtering). If scoring non-IMS data, this should be empty, otherwise a missing information exception is thrown when integrating spectra for scoring.
     * @param[out] idscores_out Reused output vectors for target or decoy identification scores
     * @param[in] mobilogram_consumer Optional consumer to write out extracted ion mobilograms
    */
    void scoreIdentification_(MRMTransitionGroupType& transition_group_identification,
                              MRMTransitionGroupType& transition_group_detection,
                              OpenSwathScoring& scorer,
                              const size_t feature_idx,
                              const Int64 feature_id,
                              const std::vector<std::string>& native_ids_detection,
                              const std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators_identification,
                              const double det_intensity_ratio_score,
                              const double det_mi_ratio_score,
                              const std::vector<OpenSwath::SwathMap>& swath_maps,
                              const double drift_target,
                              RangeMobility& im_range,
                              OpenSwath_Ind_Scores_Pooled& idscores_out,
                              MobilogramParquetConsumer* mobilogram_consumer = nullptr) const;

    void prepareScoredFeatureOutput_(OpenMS::MRMFeature& mrmfeature,
                                     const PeptideType& pep,
                                     ProteaseDigestion& pd,
                                     bool swath_present,
                                     double precursor_mz,
                                     bool ms1only) const;

    void addScoreMetaValues_(OpenMS::MRMFeature& mrmfeature,
                             const MRMTransitionGroupType& transition_group_detection,
                             const std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
                             const std::vector<OpenSwath::ISignalToNoisePtr>& ms1_signal_noise_estimators,
                             double expected_rt,
                             bool swath_present,
                             bool ms1only) const;

    void prepareFeatureOutput_(OpenMS::MRMFeature& mrmfeature, bool ms1only, int charge) const;

    /// Synchronize members with param class
    void updateMembers_() override;

    // parameters
    double rt_extraction_window_;
    double quantification_cutoff_;
    int stop_report_after_feature_;
    bool write_convex_hull_;
    bool strict_;
    bool use_ms1_ion_mobility_;
    bool apply_im_peak_picking_;
    String scoring_model_;
    String enzyme_;

    // scoring parameters
    double rt_normalization_factor_;
    int add_up_spectra_;
    String spectrum_addition_method_ ;
    String spectrum_merge_method_type_;
    double spacing_for_spectra_resampling_;
    double merge_spectra_by_peak_width_fraction_;
    double uis_threshold_sn_;
    double uis_threshold_peak_area_;

    double sn_win_len_;
    unsigned int sn_bin_count_;
    bool write_log_messages_;

    double im_extra_drift_;

    std::unordered_map<OpenMS::String, const PeptideType*> PeptideRefMap_;
    OpenSwath_Scores_Usage su_;
    OpenMS::DIAScoring diascoring_;
    OpenMS::EmgScoring emgscoring_;

    // data
    OpenSwath::SpectrumAccessPtr ms1_map_;

  };
}

#undef run_identifier
