// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
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
    for(std::map<String, int>::const_iterator it = feature_map_.begin(); it != feature_map_.end(); it++ )
    {
      result.push_back(it->first);
    }
  }
}

#endif
