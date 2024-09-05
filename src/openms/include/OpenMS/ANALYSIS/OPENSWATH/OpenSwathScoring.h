// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// data access
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

//logging
#include <OpenMS/CONCEPT/LogStream.h>

struct RangeMZ;
struct RangeMobility;

namespace OpenMS
{
  /** @brief A class that calls the scoring routines
   *
   * Use this class to invoke the individual OpenSWATH scoring routines.
   *
  */
  class OPENMS_DLLAPI OpenSwathScoring
  {
    typedef OpenSwath::LightCompound CompoundType;
    typedef OpenSwath::LightTransition TransitionType;

    enum class SpectrumAdditionMethod
    {
      ADDITION,
      RESAMPLE
    };

    double rt_normalization_factor_;
    double spacing_for_spectra_resampling_;
    int add_up_spectra_;
    SpectrumAdditionMethod spectra_addition_method_;
    double im_drift_extra_pcnt_;
    OpenSwath_Scores_Usage su_;
    bool use_ms1_ion_mobility_; ///< whether to use MS1 ion mobility extraction in DIA scores
    const std::string ION_MOBILITY_DESCRIPTION = "Ion Mobility";

  public:

    /// Constructor
    OpenSwathScoring();

    /// Destructor
    ~OpenSwathScoring();

    /** @brief Initialize the scoring object
     *
     * Sets the parameters for the scoring.
     *
     * @param rt_normalization_factor Specifies the range of the normalized retention time space
     * @param add_up_spectra How many spectra to add up (default 1)
     * @param spacing_for_spectra_resampling Spacing factor for spectra addition
     * @param drift_extra Extend the extraction window to gain a larger field of view beyond drift_upper - drift_lower (in percent)
     * @param su Which scores to actually compute
     * @param spectrum_addition_method Method to use for spectrum addition (valid: "simple", "resample")
     * @param use_ms1_ion_mobility Use MS1 ion mobility extraction in DIA scores
     *
    */
    void initialize(double rt_normalization_factor,
                    int add_up_spectra,
                    double spacing_for_spectra_resampling,
                    const double drift_extra,
                    const OpenSwath_Scores_Usage & su,
                    const std::string& spectrum_addition_method,
                    bool use_ms1_ion_mobility);

    /** @brief Score a single peakgroup in a chromatogram using only chromatographic properties.
     *
     * This function only uses the chromatographic properties (coelution,
     * signal to noise, etc.) of a peakgroup in a chromatogram to compute
     * scores. If more information is available, also consider using the
     * library based scoring and the full-spectrum based scoring.
     *
     * The scores are returned in the OpenSwath_Scores object. Only those
     * scores specified in the OpenSwath_Scores_Usage object are computed.
     *
     * @param imrmfeature The feature to be scored
     * @param native_ids The list of native ids (giving a canonical ordering of the transitions)
     * @param precursor_ids The list of precursor ids
     * @param normalized_library_intensity The weights to be used for each transition (e.g. normalized library intensities)
     * @param signal_noise_estimators The signal-to-noise estimators for each transition
     * @param scores The object to store the result
     *
    */
    void calculateChromatographicScores(OpenSwath::IMRMFeature* imrmfeature,
                                        const std::vector<std::string>& native_ids,
                                        const std::vector<std::string>& precursor_ids,
                                        const std::vector<double>& normalized_library_intensity,
                                        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
                                        OpenSwath_Scores & scores) const;

    /** @brief Score identification transitions against detection transitions of a single peakgroup
     * in a chromatogram using only chromatographic properties.
     *
     * This function only uses the chromatographic properties (coelution,
     * signal to noise, etc.) of a peakgroup in a chromatogram to compute
     * scores. The scores are computed by scoring identification against detection
     * transitions.
     *
     * The scores are returned in the OpenSwath_Scores object. Only those
     * scores specified in the OpenSwath_Scores_Usage object are computed.
     *
     * @param imrmfeature The feature to be scored
     * @param native_ids_identification The list of identification native ids (giving a canonical ordering of the transitions)
     * @param native_ids_detection The list of detection native ids (giving a canonical ordering of the transitions)
     * @param signal_noise_estimators The signal-to-noise estimators for each transition
     * @param scores The object to store the result
     *
    */
    void calculateChromatographicIdScores(OpenSwath::IMRMFeature* imrmfeature,
                                          const std::vector<std::string>& native_ids_identification,
                                          const std::vector<std::string>& native_ids_detection,
                                          std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
                                          OpenSwath_Ind_Scores & scores) const;

    /** @brief Score a single chromatographic feature against a spectral library
     *
     * The spectral library is provided in a set of transition objects and a
     * peptide object. Both contain information about the expected elution time
     * on the chromatography and the relative intensity of the transitions.
     *
     * The scores are returned in the OpenSwath_Scores object.
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param compound The compound corresponding to the library transitions
     * @param normalized_feature_rt The retention time of the feature in normalized space
     * @param scores The object to store the result
     *
    */
    void calculateLibraryScores(OpenSwath::IMRMFeature* imrmfeature,
                                const std::vector<TransitionType> & transitions,
                                const CompoundType& compound,
                                const double normalized_feature_rt,
                                OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object.
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param swath_maps The SWATH-MS (DIA) maps from which to retrieve full MS/MS spectra at the chromatographic peak apices
     * @param ms1_map The corresponding MS1 (precursor ion map) from which the precursor spectra can be retrieved (optional, may be NULL)
     * @param diascoring DIA Scoring object to use for scoring
     * @param compound The compound corresponding to the library transitions
     * @param scores The object to store the result
     * @param mzerror_ppm m/z and mass error (in ppm) for all transitions
     * @param[in] drift_target target drift value
     * @param[in] range_im drift time lower and upper bounds
     *
    */
    void calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature,
                            const std::vector<TransitionType>& transitions,
                            const std::vector<OpenSwath::SwathMap>& swath_maps,
                            const OpenSwath::SpectrumAccessPtr& ms1_map,
                            const OpenMS::DIAScoring& diascoring,
                            const CompoundType& compound,
                            OpenSwath_Scores& scores,
                            std::vector<double>& mzerror_ppm,
                            const double drift_target,
                            const RangeMobility& range_im);

    /** @brief Score a single chromatographic feature using the precursor map.
     *
     * The scores are returned in the OpenSwath_Scores object.
     *
     * @param ms1_map The MS1 (precursor ion map) from which the precursor spectra can be retrieved
     * @param diascoring DIA Scoring object to use for scoring
     * @param precursor_mz The m/z ratio of the precursor
     * @param rt The compound retention time
     * @param compound the compound sequence
     * @param im_range drift time lower and upper bounds
     * @param scores The object to store the result
     *
    */
    void calculatePrecursorDIAScores(const OpenSwath::SpectrumAccessPtr& ms1_map,
                                     const OpenMS::DIAScoring& diascoring,
                                     double precursor_mz,
                                     double rt,
                                     const CompoundType& compound,
                                     RangeMobility im_range,
                                     OpenSwath_Scores& scores);

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object.
     *
     * @param imrmfeature The feature to be scored
     * @param transition The library transition to score the feature against
     * @param swath_maps The SWATH-MS (DIA) maps from which to retrieve full MS/MS spectra at the chromatographic peak apices
     * @param range_im drift time lower and upper bounds
     * @param diascoring DIA Scoring object to use for scoring
     * @param scores The object to store the result
     *
    */
    void calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature,
                              const TransitionType & transition,
                              const std::vector<OpenSwath::SwathMap>& swath_maps,
                              RangeMobility& range_im,
                              const OpenMS::DIAScoring & diascoring,
                              OpenSwath_Scores & scores);

    /** @brief Computing the normalized library intensities from the transition objects
     *
     * The intensities are normalized such that the sum to one.
     *
     * @param[in] transitions The library transition to score the feature against
     * @param[out] normalized_library_intensity The resulting normalized library intensities
     *
    */
    void getNormalized_library_intensities_(const std::vector<TransitionType> & transitions,
                                            std::vector<double>& normalized_library_intensity);

    /** @brief Prepares a spectrum for DIA analysis (single map)
     *
     * This function will fetch a vector of spectrum pointers to be used in DIA analysis.
     * If nr_spectra_to_add == 1, then a vector of length 1 will be returned
     *
     *   - Case \#1: Non SONAR data and "simple" addition selected - Array of length "nr_spectra_to_add" returned corresponding with "nr_spectra_to_add" spectra
     *   - Case \#2: Non SONAR data and "resampling addition selected - Array of length 1 of the resampled spectrum returned
     *   - Case \#3: SONAR data - Array of length 1 containing the added/resampled spectrum returned
     *
     * For case \#2 and \#3 result is
     * all spectra summed up (add) with the intensities of multiple spectra a single
     * swath map (assuming these are regular SWATH / DIA maps) around the given
     * retention time and return an "averaged" spectrum which may contain less noise.
     *
     * For case \#1 this processing is done downstream in DIA scores to speed up computation time
     *
     * @param[in] swath_maps The map containing the spectra
     * @param[in] RT The target retention time
     * @param[in] nr_spectra_to_add How many spectra to add up
     * @param[in] im_range Drift time lower and upper bounds
     * @return Vector of spectra to be used
     *
    */
    SpectrumSequence fetchSpectrumSwath(std::vector<OpenSwath::SwathMap> swath_maps, double RT, int nr_spectra_to_add, const RangeMobility& im_range);


   /** @brief Prepares a spectrum for DIA analysis (multiple map)
     *
     * This function will fetch a SpectrumSequence to be used in DIA analysis.
     * If nr_spectra_to_add == 1, then a vector of length 1 will be returned.
     * Spectra are prepared differently based on the condition
     * Case #1: Non SONAR data and "simple" addition selected - Array of length "nr_spectra_to_add" returned corresponding with "nr_spectra_to_add" spectra
     * Case #2: Non SONAR data and "resampling addition selected - Array of length 1 of the resampled spectrum returned
     * Case #3: SONAR data - Array of length 1 containing the added/resampled spectrum returned
     *
     * For case #2 and #3 result is
     * all spectra summed up (add) with the intensities of multiple spectra a single
     * swath map (assuming these are regular SWATH / DIA maps) around the given
     * retention time and return an "averaged" spectrum which may contain less noise.
     * Spectra are also filtered and summed across drift time to transform an ion mobility spectrum into a non ion mobility spectrum
     *
     * For case #1 this processing is done downstream in DIA scores to speed up computation time, furthermore drift time filtering is done downstream (these parameters are ignored)
     *
     * @param[in] swath_map The map containing the spectra
     * @param[in] RT The target retention time
     * @param[in] nr_spectra_to_add How many spectra to add up
     * @param[in] im_range mobility range, only used if resampling spectrum addition chosen
     *
     * @return Vector of spectra to be used
     *
    */
    SpectrumSequence fetchSpectrumSwath(OpenSwath::SpectrumAccessPtr swath_map, double RT, int nr_spectra_to_add, const RangeMobility& im_range);
  };
}