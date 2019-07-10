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
    pg_scores_(rhs.pg_scores_),
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
    return pg_scores_;
  }

  OpenSwath_Scores & MRMFeature::getScores()
  {
    return pg_scores_;
  }

  void MRMFeature::setScores(const OpenSwath_Scores & scores)
  {
    pg_scores_ = scores;
  }

  void MRMFeature::addScore(const String & score_name, double score)
  {
    setMetaValue(score_name, score);
  }

  void MRMFeature::addFeature(const Feature & feature, const String& key)
  {
    features_.push_back(feature);
    feature_map_[key] = Int(features_.size()) - 1;
  }

  void MRMFeature::addFeature(Feature && feature, const String& key)
  {
    features_.push_back(std::move(feature));
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

  void MRMFeature::addPrecursorFeature(Feature && feature, const String& key)
  {
    precursor_features_.push_back(std::move(feature));
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

    setMetaValue(id + "transition_names", idscores.ind_transition_names);
    setMetaValue(id + "num_transitions", idscores.ind_num_transitions);
    setMetaValue(id + "area_intensity", idscores.ind_area_intensity);
    setMetaValue(id + "total_area_intensity", idscores.ind_total_area_intensity);
    setMetaValue(id + "intensity_score", idscores.ind_intensity_score);
    setMetaValue(id + "intensity_ratio_score", idscores.ind_intensity_ratio);
    setMetaValue(id + "apex_intensity", idscores.ind_apex_intensity);
    setMetaValue(id + "total_mi", idscores.ind_total_mi);
    setMetaValue(id + "transition_names", idscores.ind_transition_names);
    setMetaValue(id + "ind_log_intensity", idscores.ind_log_intensity);
    setMetaValue(id + "ind_xcorr_coelution", idscores.ind_xcorr_coelution_score);
    setMetaValue(id + "ind_xcorr_shape", idscores.ind_xcorr_shape_score);
    setMetaValue(id + "ind_log_sn_score", idscores.ind_log_sn_score);
    setMetaValue(id + "ind_isotope_correlation", idscores.ind_isotope_correlation);
    setMetaValue(id + "ind_isotope_overlap", idscores.ind_isotope_overlap);
    setMetaValue(id + "ind_massdev_score", idscores.ind_massdev_score);
    setMetaValue(id + "ind_mi_score", idscores.ind_mi_score);
    setMetaValue(id + "ind_mi_ratio_score", idscores.ind_mi_ratio);
  }
}

