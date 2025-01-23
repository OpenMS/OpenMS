// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <vector>

namespace OpenMS
{

  /** @brief A structure to store which scores should be used by the OpenSWATH Algorithm
   *
   * This can be used to turn on/off individual scores.
  */
  struct OPENMS_DLLAPI OpenSwath_Scores_Usage
  {
    bool use_coelution_score_ = true;
    bool use_shape_score_ = true;
    bool use_rt_score_ = true;
    bool use_library_score_ = true;
    bool use_elution_model_score_ = true;
    bool use_intensity_score_ = true;
    bool use_total_xic_score_ = true;
    bool use_total_mi_score_ = true;
    bool use_nr_peaks_score_ = true;
    bool use_sn_score_ = true;
    bool use_mi_score_ = true;
    bool use_dia_scores_ = true;
    bool use_im_scores = true;
    bool use_ms1_correlation = true;
    bool use_ms1_fullscan = true;
    bool use_ms1_mi = true;
    bool use_uis_scores = true;
    bool use_ionseries_scores = true;
    bool use_ms2_isotope_scores = true;
    bool use_peak_shape_metrics = false;
  };

  /** @brief A structure to hold the different scores computed by OpenSWATH
   *
   * This struct is used to store the individual OpenSWATH (sub-)scores. It
   * also allows to compute some preliminary quality score for a feature by
   * using a predefined combination of the individual scores determined using
   * LDA.
   *
  */
  struct OPENMS_DLLAPI OpenSwath_Scores
  {
    double elution_model_fit_score = 0;
    double library_corr = 0;
    double library_norm_manhattan = 0;
    double library_rootmeansquare = 0;
    double library_sangle = 0;
    double norm_rt_score = 0;

    double isotope_correlation = 0;
    double isotope_overlap = 0;
    double massdev_score = 0;
    double xcorr_coelution_score = 0;
    double xcorr_shape_score = 0;

    double yseries_score = 0;
    double bseries_score = 0;
    double log_sn_score = 0;

    double weighted_coelution_score = 0;
    double weighted_xcorr_shape = 0;
    double weighted_massdev_score = 0;

    double ms1_xcorr_coelution_score = -1;
    double ms1_xcorr_coelution_contrast_score = 0;
    double ms1_xcorr_coelution_combined_score = 0;
    double ms1_xcorr_shape_score = -1;
    double ms1_xcorr_shape_contrast_score = 0;
    double ms1_xcorr_shape_combined_score = 0;
    double ms1_ppm_score = 0;
    double ms1_isotope_correlation = 0;
    double ms1_isotope_overlap = 0;
    double ms1_mi_score = -1;
    double ms1_mi_contrast_score = 0;
    double ms1_mi_combined_score = 0;

    double im_xcorr_coelution_score = 0;
    double im_xcorr_shape_score = 0;
    double im_delta_score = 0;
    double im_ms1_delta_score = 0;
    double im_drift = 0;
    double im_drift_weighted = 0;
    double im_delta = -1;
    double im_ms1_contrast_coelution = 0;
    double im_ms1_contrast_shape = 0;
    double im_ms1_sum_contrast_coelution = 0;
    double im_ms1_sum_contrast_shape = 0;
    double im_ms1_drift = 0;
    double im_ms1_delta = -1;

    double library_manhattan = 0;
    double library_dotprod = 0;
    double intensity = 0;
    double total_xic = 0;
    double nr_peaks = 0;
    double sn_ratio = 0;
    double mi_score = 0;
    double weighted_mi_score = 0;

    double rt_difference = 0;
    double normalized_experimental_rt = 0;
    double raw_rt_score = 0;

    double dotprod_score_dia = 0;
    double manhatt_score_dia = 0;

    OpenSwath_Scores() = default;

    double get_quick_lda_score(double library_corr_,
                               double library_norm_manhattan_,
                               double norm_rt_score_,
                               double xcorr_coelution_score_,
                               double xcorr_shape_score_,
                               double log_sn_score_) const;

    /** @brief A quick LDA model based non-DIA scores
     *
     * A quick model LDA average model on which uses only non-DIA scores
     * (library_corr, library_norm_manhattan, norm_rt_score,
     * xcorr_coelution_score, xcorr_shape_score, log_sn_score and
     * elution_model_fit_score).
     *
     * @returns A score which is better when more negative
     *
    */
    double calculate_lda_prescore(const OpenSwath_Scores& scores) const;

    /** @brief A scoring model for peak groups with a single transition
     *
     * Manually derived scoring model for single transition peakgroups, only
     * uses norm_rt_score, log_sn_score, and elution_model_fit_score.
     *
     * @returns A score which is better when more negative
     *
    */
    double calculate_lda_single_transition(const OpenSwath_Scores& scores) const;

    /** @brief A full LDA model using DIA and non-DIA scores
     *
     * A LDA average model which uses all available scores.
     *
     * @returns A score which is better when more negative
     *
    */
    double calculate_swath_lda_prescore(const OpenSwath_Scores& scores) const;

  };

  struct OPENMS_DLLAPI OpenSwath_Ind_Scores
  {
    int ind_num_transitions = 0;
    std::vector<OpenMS::String> ind_transition_names;
    std::vector<double> ind_isotope_correlation;
    std::vector<double> ind_isotope_overlap;
    std::vector<double> ind_massdev_score;
    std::vector<double> ind_xcorr_coelution_score;
    std::vector<double> ind_xcorr_shape_score;
    std::vector<double> ind_log_sn_score;
    std::vector<double> ind_area_intensity;
    std::vector<double> ind_total_area_intensity;
    std::vector<double> ind_intensity_score;
    std::vector<double> ind_apex_intensity;
    std::vector<double> ind_apex_position;
    std::vector<double> ind_fwhm;
    std::vector<double> ind_total_mi;
    std::vector<double> ind_log_intensity;
    std::vector<double> ind_intensity_ratio;
    std::vector<double> ind_mi_ratio;
    std::vector<double> ind_mi_score;
    
    // peak shape metrics
    std::vector<double> ind_start_position_at_5;
    std::vector<double> ind_end_position_at_5;
    std::vector<double> ind_start_position_at_10;
    std::vector<double> ind_end_position_at_10;
    std::vector<double> ind_start_position_at_50;
    std::vector<double> ind_end_position_at_50;
    std::vector<double> ind_total_width;
    std::vector<double> ind_tailing_factor;
    std::vector<double> ind_asymmetry_factor;
    std::vector<double> ind_slope_of_baseline;
    std::vector<double> ind_baseline_delta_2_height;
    std::vector<double> ind_points_across_baseline;
    std::vector<double> ind_points_across_half_height;

    OpenSwath_Ind_Scores() = default;

  };

}


