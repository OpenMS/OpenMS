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
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace OpenMS
{

  /** @brief A structure to store which scores should be used by the Algorithm
   *
   * This can be used to turn on/off individual scores.
  */
  struct OPENMS_DLLAPI OpenSwath_Scores_Usage
  {
    // Which scores to use
    bool use_coelution_score_;
    bool use_shape_score_;
    bool use_rt_score_;
    bool use_library_score_;
    bool use_elution_model_score_;
    bool use_intensity_score_;
    bool use_total_xic_score_;
    bool use_total_mi_score_;
    bool use_nr_peaks_score_;
    bool use_sn_score_;
    bool use_mi_score_;
    bool use_dia_scores_;
    bool use_sonar_scores;
    bool use_ms1_correlation;
    bool use_ms1_fullscan;
    bool use_ms1_mi;
    bool use_uis_scores;
    
    OpenSwath_Scores_Usage() :
      use_coelution_score_(true),
      use_shape_score_(true),
      use_rt_score_(true),
      use_library_score_(true),
      use_elution_model_score_(true),
      use_intensity_score_(true),
      use_total_xic_score_(true),
      use_total_mi_score_(true),
      use_nr_peaks_score_(true),
      use_sn_score_(true),
      use_mi_score_(true),
      use_dia_scores_(true),
      use_sonar_scores(true),
      use_ms1_correlation(true),
      use_ms1_fullscan(true),
      use_ms1_mi(true),
      use_uis_scores(true)
    {}

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
    double elution_model_fit_score;
    double library_corr;
    double library_norm_manhattan;
    double library_rootmeansquare;
    double library_sangle;
    double norm_rt_score;
    double isotope_correlation;
    std::string ind_isotope_correlation;
    double isotope_overlap;
    std::string ind_isotope_overlap;
    double massdev_score;
    std::string ind_massdev_score;
    double xcorr_coelution_score;
    std::string ind_xcorr_coelution_score;
    double xcorr_shape_score;
    std::string ind_xcorr_shape_score;
    double yseries_score;
    double bseries_score;
    double log_sn_score;
    std::string ind_log_sn_score;
    int ind_num_transitions;
    std::string ind_transition_names;
    std::string ind_area_intensity;
    std::string ind_total_area_intensity;
    std::string ind_intensity_score;
    std::string ind_apex_intensity;
    std::string ind_total_mi;
    std::string ind_log_intensity;
    std::string ind_intensity_ratio;
    std::string ind_mi_ratio;

    double weighted_coelution_score;
    double weighted_xcorr_shape;
    double weighted_massdev_score;
   
    double ms1_xcorr_coelution_score;
    double ms1_xcorr_coelution_contrast_score;
    double ms1_xcorr_coelution_combined_score;
    double ms1_xcorr_shape_score;
    double ms1_xcorr_shape_contrast_score;
    double ms1_xcorr_shape_combined_score;
    double ms1_ppm_score;
    double ms1_isotope_correlation;
    double ms1_isotope_overlap;
    double ms1_mi_score;
    double ms1_mi_contrast_score;
    double ms1_mi_combined_score;

    double sonar_sn;
    double sonar_diff;
    double sonar_trend;
    double sonar_rsq;
    double sonar_shape;
    double sonar_lag;

    double library_manhattan;
    double library_dotprod;
    double intensity;
    double total_xic;
    double nr_peaks;
    double sn_ratio;
    double mi_score;
    std::string ind_mi_score;
    double weighted_mi_score;

    double rt_difference;
    double normalized_experimental_rt;
    double raw_rt_score;

    double dotprod_score_dia;
    double manhatt_score_dia;

    OpenSwath_Scores() :
      elution_model_fit_score(0),
      library_corr(0),
      library_norm_manhattan(0),
      library_rootmeansquare(0),
      library_sangle(0),
      norm_rt_score(0),
      isotope_correlation(0),
      ind_isotope_correlation(""),
      isotope_overlap(0),
      ind_isotope_overlap(""),
      massdev_score(0),
      ind_massdev_score(""),
      xcorr_coelution_score(0),
      ind_xcorr_coelution_score(""),
      xcorr_shape_score(0),
      ind_xcorr_shape_score(""),
      yseries_score(0),
      bseries_score(0),
      log_sn_score(0),
      ind_log_sn_score(""),
      ind_num_transitions(0),
      ind_transition_names(""),
      ind_area_intensity(""),
      ind_total_area_intensity(""),
      ind_intensity_score(""),
      ind_apex_intensity(""),
      ind_total_mi(""),
      ind_log_intensity(""),
      ind_intensity_ratio(""),
      ind_mi_ratio(""),
      weighted_coelution_score(0),
      weighted_xcorr_shape(0),
      weighted_massdev_score(0),
      ms1_xcorr_coelution_score(-1),
      ms1_xcorr_coelution_contrast_score(0),
      ms1_xcorr_coelution_combined_score(0),
      ms1_xcorr_shape_score(-1),
      ms1_xcorr_shape_contrast_score(0),
      ms1_xcorr_shape_combined_score(0),
      ms1_ppm_score(0),
      ms1_isotope_correlation(0),
      ms1_isotope_overlap(0),
      ms1_mi_score(-1),
      ms1_mi_contrast_score(0),
      ms1_mi_combined_score(0),
      sonar_sn(0),
      sonar_diff(0),
      sonar_trend(0),
      sonar_rsq(0),
      sonar_shape(0),
      sonar_lag(0),
      library_manhattan(0),
      library_dotprod(0),
      intensity(0),
      total_xic(0),
      nr_peaks(0),
      sn_ratio(0),
      mi_score(0),
      ind_mi_score(""),
      weighted_mi_score(0),
      dotprod_score_dia(0),
      manhatt_score_dia(0)
    {
    }


    double get_quick_lda_score(double library_corr_, double library_norm_manhattan_, double norm_rt_score_, double xcorr_coelution_score_,
                               double xcorr_shape_score_, double log_sn_score_) const;

    double calculate_lda_prescore(const OpenSwath_Scores& scores) const;

    double calculate_lda_single_transition(const OpenSwath_Scores& scores) const;

    double calculate_swath_lda_prescore(const OpenSwath_Scores& scores) const;

  };

}


