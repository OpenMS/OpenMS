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

#include <OpenMS/KERNEL/MRMFeature.h>

namespace OpenMS
{

  MRMFeature::MRMFeature() :
    Feature()
  {
  }

  /// Copy constructor
  MRMFeature::MRMFeature(const MRMFeature & rhs) :
    Feature(rhs),
    features_(rhs.features_),
    precursor_features_(rhs.precursor_features_),
    feature_map_(rhs.feature_map_),
    precursor_feature_map_(rhs.precursor_feature_map_)
  {
    setScores(rhs.getScores());
  }

  /// Assignment operator
  MRMFeature & MRMFeature::operator = (const MRMFeature &rhs)
  {
    if (&rhs == this)
      return *this;

    Feature::operator = (rhs);
    setScores(rhs.getScores());
    features_ = rhs.features_;
    precursor_features_ = rhs.precursor_features_;
    feature_map_ = rhs.feature_map_;
    precursor_feature_map_ = rhs.precursor_feature_map_;

    return *this;
  }

  MRMFeature::~MRMFeature()
  {
  }

  const OpenSwath_Scores & MRMFeature::getScores() const
  {
    return scores_;
  }

  OpenSwath_Scores & MRMFeature::getScores()
  {
    return scores_;
  }

  void MRMFeature::setScores(const OpenSwath_Scores & scores)
  {
    scores_ = scores;
  }

  void MRMFeature::addFeature(const Feature & feature, const String& key)
  {
    features_.push_back(feature);
    feature_map_[key] = Int(features_.size()) - 1;
  }

  Feature & MRMFeature::getFeature(const String& key) 
  {
    return features_.at(feature_map_[key]);
  }

  const Feature & MRMFeature::getFeature(const String& key) const 
  {
    return features_.at(feature_map_.at(key));
  }

  const std::vector<Feature> & MRMFeature::getFeatures() const
  {
    return features_;
  }

  void MRMFeature::getFeatureIDs(std::vector<String> & result) const
  {
    for (std::map<String, int>::const_iterator it = feature_map_.begin(); it != feature_map_.end(); ++it)
    {
      result.push_back(it->first);
    }
  }

  void MRMFeature::addPrecursorFeature(const Feature & feature, const String& key)
  {
    precursor_features_.push_back(feature);
    precursor_feature_map_[key] = Int(precursor_features_.size()) - 1;
  }

  void MRMFeature::getPrecursorFeatureIDs(std::vector<String> & result) const
  {
    for (std::map<String, int>::const_iterator it = precursor_feature_map_.begin(); it != precursor_feature_map_.end(); ++it)
    {
      result.push_back(it->first);
    }
  }

  Feature & MRMFeature::getPrecursorFeature(String key)
  {
    return precursor_features_.at(precursor_feature_map_[key]);
  }

  const Feature & MRMFeature::getPrecursorFeature(String key) const
  {
    return precursor_features_.at(precursor_feature_map_.at(key));
  }

  void MRMFeature::IDScoresAsMetaValue(bool decoy, const OpenSwath_Ind_Scores& idscores)
  {
    String id = "id_target_";
    if (decoy) id = "id_decoy_";

    this->setMetaValue(id + "transition_names", idscores.ind_transition_names);
    this->setMetaValue(id + "num_transitions", idscores.ind_num_transitions);
    this->setMetaValue(id + "area_intensity", idscores.ind_area_intensity);
    this->setMetaValue(id + "total_area_intensity", idscores.ind_total_area_intensity);
    this->setMetaValue(id + "intensity_score", idscores.ind_intensity_score);
    this->setMetaValue(id + "intensity_ratio_score", idscores.ind_intensity_ratio);
    this->setMetaValue(id + "apex_intensity", idscores.ind_apex_intensity);
    this->setMetaValue(id + "total_mi", idscores.ind_total_mi);
    this->setMetaValue(id + "transition_names", idscores.ind_transition_names);
    this->setMetaValue(id + "ind_log_intensity", idscores.ind_log_intensity);
    this->setMetaValue(id + "ind_xcorr_coelution", idscores.ind_xcorr_coelution_score);
    this->setMetaValue(id + "ind_xcorr_shape", idscores.ind_xcorr_shape_score);
    this->setMetaValue(id + "ind_log_sn_score", idscores.ind_log_sn_score);
    this->setMetaValue(id + "ind_isotope_correlation", idscores.ind_isotope_correlation);
    this->setMetaValue(id + "ind_isotope_overlap", idscores.ind_isotope_overlap);
    this->setMetaValue(id + "ind_massdev_score", idscores.ind_massdev_score);
    this->setMetaValue(id + "ind_mi_score", idscores.ind_mi_score);
    this->setMetaValue(id + "ind_mi_ratio_score", idscores.ind_mi_ratio);
  }

  void MRMFeature::scoresAsMetaValue(bool ms1only, const OpenSwath_Scores_Usage& su_)
  {
      if (ms1only)
      {
        if (su_.use_sn_score_) 
        { 
          this->setMetaValue("sn_ratio", scores_.sn_ratio);
          this->setMetaValue("var_log_sn_score", scores_.log_sn_score); 
        }

        if (su_.use_rt_score_)
        {
          this->setMetaValue("delta_rt", this->getRT() - expected_rt_);
          this->setMetaValue("assay_rt", expected_rt_);
          this->setMetaValue("norm_RT", scores_.normalized_experimental_rt);
          this->setMetaValue("rt_score", scores_.raw_rt_score);
          this->setMetaValue("var_norm_rt_score", scores_.norm_rt_score);
        }

        if (su_.use_ms1_fullscan)
        {
          this->setMetaValue("var_ms1_ppm_diff", scores_.ms1_ppm_score);
          this->setMetaValue("var_ms1_isotope_correlation", scores_.ms1_isotope_correlation);
          this->setMetaValue("var_ms1_isotope_overlap", scores_.ms1_isotope_overlap);
        }
      }
      else
      {
        if (su_.use_coelution_score_)
        {
          this->setMetaValue("var_xcorr_coelution", scores_.xcorr_coelution_score);
          this->setMetaValue("var_xcorr_coelution_weighted", scores_.weighted_coelution_score);
        }
        if (su_.use_shape_score_)
        {
          this->setMetaValue("var_xcorr_shape", scores_.xcorr_shape_score);
          this->setMetaValue("var_xcorr_shape_weighted", scores_.weighted_xcorr_shape);
        }
        if (su_.use_library_score_)
        {
          this->setMetaValue("var_library_corr", scores_.library_corr);
          this->setMetaValue("var_library_rmsd", scores_.library_norm_manhattan);
          this->setMetaValue("var_library_sangle", scores_.library_sangle);
          this->setMetaValue("var_library_rootmeansquare", scores_.library_rootmeansquare);
          this->setMetaValue("var_library_manhattan", scores_.library_manhattan);
          this->setMetaValue("var_library_dotprod", scores_.library_dotprod);
        }
        if (su_.use_rt_score_)
        {
          this->setMetaValue("delta_rt", this->getRT() - expected_rt_);
          this->setMetaValue("assay_rt", expected_rt_);
          this->setMetaValue("norm_RT", scores_.normalized_experimental_rt);
          this->setMetaValue("rt_score", scores_.raw_rt_score);
          this->setMetaValue("var_norm_rt_score", scores_.norm_rt_score);
        }

        // TODO do we really want these intensity scores_ ?
        if (su_.use_intensity_score_)
        {
          if ((double)this->getMetaValue("total_xic") > 0)
          {
            this->setMetaValue("var_intensity_score", this->getIntensity() / (double)this->getMetaValue("total_xic"));
          }
          else
          {
            this->setMetaValue("var_intensity_score", 0);
          }
        }

        if (su_.use_total_xic_score_) { this->setMetaValue("total_xic", (double)this->getMetaValue("total_xic")); }
        if (su_.use_total_mi_score_) { this->setMetaValue("total_mi", (double)this->getMetaValue("total_mi")); }
        if (su_.use_nr_peaks_score_) { this->setMetaValue("nr_peaks", scores_.nr_peaks); }

        if (su_.use_sn_score_)
        {
          this->setMetaValue("sn_ratio", scores_.sn_ratio);
          this->setMetaValue("var_log_sn_score", scores_.log_sn_score);
        }

        if (su_.use_mi_score_)
        {
          this->setMetaValue("var_mi_score", scores_.mi_score);
          this->setMetaValue("var_mi_weighted_score", scores_.weighted_mi_score);
          if (su_.use_total_mi_score_)
          {
            if (((double)this->getMetaValue("total_mi")) > 0)
            {
              this->setMetaValue("var_mi_ratio_score", scores_.mi_score  / (double)this->getMetaValue("total_mi"));
            }
            else
            {
              this->setMetaValue("var_mi_ratio_score", 0);
            }
          }
        }

        // TODO get it working with imrmfeature
        if (su_.use_elution_model_score_)
        {
          this->setMetaValue("var_elution_model_fit_score", scores_.elution_model_fit_score);
        }

        // Add the DIA / SWATH scores_
        if (su_.use_dia_scores_)
        {
          this->setMetaValue("var_isotope_correlation_score", scores_.isotope_correlation);
          this->setMetaValue("var_isotope_overlap_score", scores_.isotope_overlap);
          this->setMetaValue("var_massdev_score", scores_.massdev_score);
          this->setMetaValue("var_massdev_score_weighted", scores_.weighted_massdev_score);
          this->setMetaValue("var_bseries_score", scores_.bseries_score);
          this->setMetaValue("var_yseries_score", scores_.yseries_score);
          this->setMetaValue("var_dotprod_score", scores_.dotprod_score_dia);
          this->setMetaValue("var_manhatt_score", scores_.manhatt_score_dia);
          if (su_.use_ms1_correlation)
          {
            if (scores_.ms1_xcorr_shape_score > -1)
            {
              this->setMetaValue("var_ms1_xcorr_shape", scores_.ms1_xcorr_shape_score);
            }
            if (scores_.ms1_xcorr_coelution_score > -1)
            {
              this->setMetaValue("var_ms1_xcorr_coelution", scores_.ms1_xcorr_coelution_score);
            }
            this->setMetaValue("var_ms1_xcorr_shape_contrast", scores_.ms1_xcorr_shape_contrast_score);
            this->setMetaValue("var_ms1_xcorr_shape_combined", scores_.ms1_xcorr_shape_combined_score);
            this->setMetaValue("var_ms1_xcorr_coelution_contrast", scores_.ms1_xcorr_coelution_contrast_score);
            this->setMetaValue("var_ms1_xcorr_coelution_combined", scores_.ms1_xcorr_coelution_combined_score);
          }
          if (su_.use_ms1_mi)
          {
            if (scores_.ms1_mi_score > -1)
            {
              this->setMetaValue("var_ms1_mi_score", scores_.ms1_mi_score);
            }
            this->setMetaValue("var_ms1_mi_contrast_score", scores_.ms1_mi_contrast_score);
            this->setMetaValue("var_ms1_mi_combined_score", scores_.ms1_mi_combined_score);
          }
          if (su_.use_ms1_fullscan)
          {
            this->setMetaValue("var_ms1_ppm_diff", scores_.ms1_ppm_score);
            this->setMetaValue("var_ms1_isotope_correlation", scores_.ms1_isotope_correlation);
            this->setMetaValue("var_ms1_isotope_overlap", scores_.ms1_isotope_overlap);
          }
        }

        if (su_.use_sonar_scores)
        {
          // set all scores_ less than 1 to zero (do not over-punish large negative scores_)
          double log_sn = 0;
          if (scores_.sonar_sn > 1) log_sn = std::log(scores_.sonar_sn);
          double log_trend = 0;
          if (scores_.sonar_trend > 1) log_trend = std::log(scores_.sonar_trend);
          double log_diff = 0;
          if (scores_.sonar_diff > 1) log_diff = std::log(scores_.sonar_diff);

          this->setMetaValue("var_sonar_lag", scores_.sonar_lag);
          this->setMetaValue("var_sonar_shape", scores_.sonar_shape);
          this->setMetaValue("var_sonar_log_sn", log_sn);
          this->setMetaValue("var_sonar_log_diff", log_diff);
          this->setMetaValue("var_sonar_log_trend", log_trend);
          this->setMetaValue("var_sonar_rsq", scores_.sonar_rsq);
        }
      }
  }

}

