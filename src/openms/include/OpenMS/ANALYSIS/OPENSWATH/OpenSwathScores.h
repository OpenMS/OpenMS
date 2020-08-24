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
    bool use_sonar_scores = true;
    bool use_ms1_correlation = true;
    bool use_ms1_fullscan = true;
    bool use_ms1_mi = true;
    bool use_uis_scores = true;
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

    double sonar_sn = 0;
    double sonar_diff = 0;
    double sonar_trend = 0;
    double sonar_rsq = 0;
    double sonar_shape = 0;
    double sonar_lag = 0;

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
    std::vector<double> ind_total_mi;
    std::vector<double> ind_log_intensity;
    std::vector<double> ind_intensity_ratio;
    std::vector<double> ind_mi_ratio;
    std::vector<double> ind_mi_score;

    OpenSwath_Ind_Scores() = default;

  };

}


