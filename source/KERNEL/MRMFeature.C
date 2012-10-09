// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_KERNEL_MRMFEATURE_C
#define OPENMS_KERNEL_MRMFEATURE_C

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
    pg_scores_(rhs.pg_scores_),
    feature_map_(rhs.feature_map_)
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
    feature_map_ = rhs.feature_map_;
    features_ = rhs.features_;

    return *this;
  }

  MRMFeature::~MRMFeature()
  {
  }

  const MRMFeature::PGScoresType & MRMFeature::getScores() const
  {
    return pg_scores_;
  }

  double MRMFeature::getScore(const String & score_name)
  {
    return pg_scores_[score_name];
  }

  void MRMFeature::setScores(const PGScoresType & scores)
  {

    for (MRMFeature::PGScoresType::const_iterator score = scores.begin();
         score != scores.end(); score++)
    {
      addScore(score->first, score->second);
    }
  }

  void MRMFeature::addScore(const String & score_name, double score)
  {
    pg_scores_[score_name] = score;
    setMetaValue(score_name, score);
  }

  void MRMFeature::addFeature(Feature & feature, String key)
  {
    features_.push_back(feature);
    feature_map_[key] = features_.size() - 1;
  }

  Feature & MRMFeature::getFeature(String key) 
  {
    return features_.at(feature_map_[key]);
  }

  const std::vector<Feature> & MRMFeature::getFeatures() const
  {
    return features_;
  }

  void MRMFeature::getFeatureIDs(std::vector<String> & result) const
  {
    for (std::map<String, int>::const_iterator it = feature_map_.begin(); it != feature_map_.end(); it++ )
    {
      result.push_back(it->first);
    }
  }
}

#endif
