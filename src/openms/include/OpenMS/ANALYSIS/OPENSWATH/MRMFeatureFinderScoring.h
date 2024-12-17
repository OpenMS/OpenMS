// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#ifdef _OPENMP
#include <omp.h>
#endif

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left, const std::pair<double, double>& right);

namespace OpenMS
{

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
     * @param chromatograms The input chromatograms
     * @param output The output features with corresponding scores
     * @param transition_exp The transition list describing the experiment
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the transition list
     * @param swath_map Optional SWATH-MS (DIA) map corresponding from which the chromatograms were extracted
     *
    */
    void pickExperiment(const PeakMap & chromatograms,
                        FeatureMap& output,
                        const TargetedExperiment& transition_exp,
                        const TransformationDescription& trafo,
                        const PeakMap& swath_map);

    /** @brief Pick and score features in a single experiment from chromatograms
     *
     * @param input The input chromatograms
     * @param output The output features with corresponding scores
     * @param transition_exp The transition list describing the experiment
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param swath_maps Optional SWATH-MS (DIA) map corresponding from which
     *                  the chromatograms were extracted. Use empty map if no
     *                  data is available.
     * @param transition_group_map Output mapping of transition groups
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
     * @param transition_exp The transition list describing the experiment
     *
    */
    void prepareProteinPeptideMaps_(const OpenSwath::LightTargetedExperiment& transition_exp);
    //@}

    /** @brief Score all peak groups of a transition group
     *
     * Iterate through all features found along the chromatograms of the
     * transition group and score each one individually.
     *
     * @param transition_group The MRMTransitionGroup to be scored (input)
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param swath_maps Optional SWATH-MS (DIA) map corresponding from which
     *                   the chromatograms were extracted. Use empty map if no
     *                   data is available.
     * @param output The output features with corresponding scores (the found
     *               features will be added to this FeatureMap).
     * @param ms1only Whether to only do MS1 scoring and skip all MS2 scoring
     *
    */
    void scorePeakgroups(MRMTransitionGroupType& transition_group,
                         const TransformationDescription & trafo,
                         const std::vector<OpenSwath::SwathMap>& swath_maps,
                         FeatureMap& output,
                         bool ms1only = false) const;

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
     * @param ms1_map The raw mass spectrometric MS1 data
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
     * @param input The input chromatograms
     * @param transition_exp The transition list describing the experiment
     * @param transition_group_map Mapping of transition groups
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param rt_extraction_window The used retention time extraction window
     *
    */
    void mapExperimentToTransitionList(const OpenSwath::SpectrumAccessPtr& input,
                                       const OpenSwath::LightTargetedExperiment& transition_exp,
                                       TransitionGroupMapType& transition_group_map,
                                       TransformationDescription trafo,
                                       double rt_extraction_window);
private:

    /** @brief Splits combined transition groups into detection transition groups
     *
     * For standard assays, transition_group_detection is identical to transition_group and the others are empty.
     *
     * @param transition_group Containing all detecting, identifying transitions
     * @param transition_group_detection To be filled with detecting transitions
    */
    void splitTransitionGroupsDetection_(const MRMTransitionGroupType& transition_group,
                                         MRMTransitionGroupType& transition_group_detection) const;

    /** @brief Splits combined transition groups into identification transition groups
     *
     * For standard assays, transition_group_identification is empty. When UIS scoring
     * is enabled, it contains the corresponding identification transitions.
     *
     * @param transition_group Containing all detecting, identifying transitions
     * @param transition_group_identification To be filled with identifying transitions
     * @param transition_group_identification_decoy To be filled with identifying decoy transitions
    */
    void splitTransitionGroupsIdentification_(const MRMTransitionGroupType& transition_group,
                                              MRMTransitionGroupType& transition_group_identification,
                                              MRMTransitionGroupType& transition_group_identification_decoy) const;

    /** @brief Provides scoring for target and decoy identification against detecting transitions
     *
     * The function is used twice, for target and decoy identification transitions. The results are
     * reported analogously to the ones for detecting transitions but must be stored separately.
     *
     * @param transition_group_identification Containing all detecting and identifying transitions
     * @param scorer An instance of OpenSwathScoring
     * @param feature_idx The index of the current feature
     * @param native_ids_detection The native IDs of the detecting transitions
     * @param det_intensity_ratio_score The intensity score of the detection transitions for normalization
     * @param det_mi_ratio_score The MI score of the detection transitions for normalization
     * @param swath_maps Optional SWATH-MS (DIA) map corresponding from which
     *                  the chromatograms were extracted. Use empty map if no
     *                  data is available.
     * @return a struct of type OpenSwath_Ind_Scores containing either target or decoy values
    */
    OpenSwath_Ind_Scores scoreIdentification_(MRMTransitionGroupType& transition_group_identification,
                                              OpenSwathScoring& scorer,
                                              const size_t feature_idx,
                                              const std::vector<std::string> & native_ids_detection,
                                              const double det_intensity_ratio_score,
                                              const double det_mi_ratio_score,
                                              const std::vector<OpenSwath::SwathMap>& swath_maps) const;

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
    String scoring_model_;

    // scoring parameters
    double rt_normalization_factor_;
    int add_up_spectra_;
    String spectrum_addition_method_ ;
    double spacing_for_spectra_resampling_;
    double uis_threshold_sn_;
    double uis_threshold_peak_area_;

    double sn_win_len_;
    unsigned int sn_bin_count_;
    bool write_log_messages_;

    double im_extra_drift_;

    // members
    std::map<OpenMS::String, const PeptideType*> PeptideRefMap_;
    OpenSwath_Scores_Usage su_;
    OpenMS::DIAScoring diascoring_;
    OpenMS::EmgScoring emgscoring_;

    // data
    OpenSwath::SpectrumAccessPtr ms1_map_;

  };
}

#undef run_identifier

