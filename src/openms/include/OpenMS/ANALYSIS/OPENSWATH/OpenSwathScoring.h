// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

    double rt_normalization_factor_;
    double spacing_for_spectra_resampling_;
    int add_up_spectra_;
    std::string spectra_addition_method_;
    double im_drift_extra_pcnt_;
    OpenSwath_Scores_Usage su_;

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
     * @param su Which scores to actually compute
     * @param spectrum_addition_method Method to use for spectrum addition (valid: "simple", "resample")
     *
    */
    void initialize(double rt_normalization_factor,
                    int add_up_spectra,
                    double spacing_for_spectra_resampling,
                    const double drift_extra,
                    const OpenSwath_Scores_Usage & su,
                    const std::string& spectrum_addition_method);

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
                                        OpenSwath_Scores & scores);

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
                                          OpenSwath_Ind_Scores & scores);

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
     * @param pep The peptide corresponding to the library transitions
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
     * @param pep The peptide corresponding to the library transitions
     * @param scores The object to store the result
     * @param mzerror_ppm m/z and mass error (in ppm) for all transitions
     * @param drift_lower Drift time lower extraction boundary
     * @param drift_upper Drift time upper extraction boundary
     *
    */
    void calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature,
                            const std::vector<TransitionType>& transitions,
                            const std::vector<OpenSwath::SwathMap>& swath_maps,
                            OpenSwath::SpectrumAccessPtr ms1_map,
                            OpenMS::DIAScoring& diascoring,
                            const CompoundType& compound,
                            OpenSwath_Scores& scores,
                            std::vector<double>& mzerror_ppm,
                            const double drift_lower,
                            const double drift_upper,
                            const double drift_target);

    /** @brief Score a single chromatographic feature using the precursor map.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param ms1_map The MS1 (precursor ion map) from which the precursor spectra can be retrieved
     * @param diascoring DIA Scoring object to use for scoring
     * @param precursor_mz The m/z ratio of the precursor
     * @param rt The compound retention time
     * @param scores The object to store the result
     * @param drift_lower Drift time lower extraction boundary
     * @param drift_upper Drift time upper extraction boundary
     *
    */
    void calculatePrecursorDIAScores(OpenSwath::SpectrumAccessPtr ms1_map, 
                                     OpenMS::DIAScoring& diascoring, 
                                     double precursor_mz, 
                                     double rt, 
                                     const CompoundType& compound, 
                                     OpenSwath_Scores& scores,
                                     double drift_lower,
                                     double drift_upper);

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param swath_maps The SWATH-MS (DIA) maps from which to retrieve full MS/MS spectra at the chromatographic peak apices
     * @param diascoring DIA Scoring object to use for scoring
     * @param scores The object to store the result
     * @param drift_lower Drift time lower extraction boundary
     * @param drift_upper Drift time upper extraction boundary
     *
    */
    void calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature,
                              const TransitionType & transition,
                              const std::vector<OpenSwath::SwathMap> swath_maps,
                              OpenMS::DIAScoring & diascoring,
                              OpenSwath_Scores & scores,
                              double drift_lower,
                              double drift_upper);

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

    /** @brief Prepares a spectrum for DIA analysis (multiple map)
     *
     * This function will sum up (add) the intensities of multiple spectra from
     * multiple swath maps (assuming these are SONAR maps of shifted precursor
     * isolation windows) around the given retention time and return an
     * "averaged" spectrum which may contain less noise.
     *
     * @param[in] swath_maps The map(s) containing the spectra
     * @param[in] RT The target retention time
     * @param[in] nr_spectra_to_add How many spectra to add up
     * @param drift_lower Drift time lower extraction boundary
     * @param drift_upper Drift time upper extraction boundary
     *
     * @return Added up spectrum
     *
    */
    OpenSwath::SpectrumPtr fetchSpectrumSwath(std::vector<OpenSwath::SwathMap> swath_maps,
                                              double RT,
                                              int nr_spectra_to_add,
                                              const double drift_lower,
                                              const double drift_upper);
    
    /** @brief Prepares a spectrum for DIA analysis (single map)
     *
     * This function will sum up (add) the intensities of multiple spectra a single
     * swath map (assuming these are regular SWATH / DIA maps) around the given 
     * retention time and return an "averaged" spectrum which may contain less noise.
     *
     * @param[in] swath_map The map containing the spectra
     * @param[in] RT The target retention time
     * @param[in] nr_spectra_to_add How many spectra to add up
     * @param drift_lower Drift time lower extraction boundary
     * @param drift_upper Drift time upper extraction boundary
     *
     * @return Added up spectrum
     *
    */
    OpenSwath::SpectrumPtr fetchSpectrumSwath(OpenSwath::SpectrumAccessPtr swath_map,
                                              double RT,
                                              int nr_spectra_to_add,
                                              const double drift_lower,
                                              const double drift_upper);

  protected:

    /** @brief Returns an averaged spectrum
     *
     * This function will sum up (add) the intensities of multiple spectra
     * around the given retention time and return an "averaged" spectrum which
     * may contain less noise.
     *
     * @param[in] swath_map The map containing the spectra
     * @param[in] RT The target retention time
     * @param[in] nr_spectra_to_add How many spectra to add up
     * @param drift_lower Drift time lower extraction boundary
     * @param drift_upper Drift time upper extraction boundary
     *
     * @return Added up spectrum
    */
    OpenSwath::SpectrumPtr getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map,
                                            double RT,
                                            int nr_spectra_to_add,
                                            const double drift_lower,
                                            const double drift_upper);

  };
}

