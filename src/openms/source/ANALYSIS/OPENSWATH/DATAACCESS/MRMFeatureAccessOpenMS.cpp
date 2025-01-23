// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

namespace OpenMS
{

  MRMFeatureOpenMS::MRMFeatureOpenMS(MRMFeature& mrmfeature) :
    mrmfeature_(mrmfeature)
  {
    std::vector<String> ids;
    mrmfeature.getFeatureIDs(ids);
    for (std::vector<String>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      boost::shared_ptr<FeatureOpenMS> ptr = boost::shared_ptr<FeatureOpenMS>(new FeatureOpenMS(mrmfeature.getFeature(*it)));
      features_[*it] = ptr;
    }

    std::vector<String> p_ids;
    mrmfeature.getPrecursorFeatureIDs(p_ids);
    for (std::vector<String>::iterator it = p_ids.begin(); it != p_ids.end(); ++it)
    {
      boost::shared_ptr<FeatureOpenMS> ptr = boost::shared_ptr<FeatureOpenMS>(new FeatureOpenMS(mrmfeature.getPrecursorFeature(*it)));
      precursor_features_[*it] = ptr;
    }
  }

  FeatureOpenMS::FeatureOpenMS(Feature& feature)
  {
    feature_ = &feature; // store raw ptr to the feature
  }

  FeatureOpenMS::~FeatureOpenMS() = default;

  void FeatureOpenMS::getRT(std::vector<double>& rt) const
  {
    OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "There needs to exactly one convex hull per feature.");
    ConvexHull2D::PointArrayType data_points = feature_->getConvexHulls()[0].getHullPoints();
    for (ConvexHull2D::PointArrayType::iterator it = data_points.begin(); it != data_points.end(); ++it)
    {
      rt.push_back(it->getX());
    }
  }

  void FeatureOpenMS::getIntensity(std::vector<double>& intens) const
  {
    OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "There needs to exactly one convex hull per feature.");
    ConvexHull2D::PointArrayType data_points = feature_->getConvexHulls()[0].getHullPoints();
    for (ConvexHull2D::PointArrayType::iterator it = data_points.begin(); it != data_points.end(); ++it)
    {
      intens.push_back(it->getY());
    }
  }

  float FeatureOpenMS::getIntensity() const
  {
    return feature_->getIntensity();
  }

  double FeatureOpenMS::getRT() const
  {
    return feature_->getRT();
  }

  MRMFeatureOpenMS::~MRMFeatureOpenMS() = default;

  boost::shared_ptr<OpenSwath::IFeature> MRMFeatureOpenMS::getFeature(std::string nativeID)
  {
    OPENMS_PRECONDITION(features_.find(nativeID) != features_.end(), "Feature needs to exist");
    return boost::static_pointer_cast<OpenSwath::IFeature>(features_[nativeID]);
  }

  boost::shared_ptr<OpenSwath::IFeature> MRMFeatureOpenMS::getPrecursorFeature(std::string nativeID)
  {
    OPENMS_PRECONDITION(precursor_features_.find(nativeID) != precursor_features_.end(), "Precursor feature needs to exist");
    return boost::static_pointer_cast<OpenSwath::IFeature>(precursor_features_[nativeID]);
  }

  std::vector<std::string> MRMFeatureOpenMS::getNativeIDs() const
  {
    std::vector<std::string> v;
    for (std::map<std::string, boost::shared_ptr<FeatureOpenMS> >::const_iterator it = features_.begin(); it != features_.end(); ++it)
    {
      v.push_back(it->first);
    }
    return v;
  }

  std::vector<std::string> MRMFeatureOpenMS::getPrecursorIDs() const
  {
    std::vector<std::string> v;
    for (std::map<std::string, boost::shared_ptr<FeatureOpenMS> >::const_iterator it = precursor_features_.begin(); it != precursor_features_.end(); ++it) 
    {
      v.push_back(it->first);
    }
    return v;
  }

  float MRMFeatureOpenMS::getIntensity() const
  {
    return mrmfeature_.getIntensity();
  }

  double MRMFeatureOpenMS::getRT() const
  {
    return mrmfeature_.getRT();
  }

  size_t MRMFeatureOpenMS::size() const
  {
    return features_.size();
  }

  // default instances
  MSSpectrum spec;
  MSChromatogram chrom;
  SignalToNoiseOpenMS< MSSpectrum> spec_signal_to_noise_openms(spec, 1.0, 3, true);
  SignalToNoiseOpenMS< MSChromatogram > chrom_signal_to_noise_openms(chrom, 1.0, 3, true);

  MRMTransitionGroup<MSSpectrum, ReactionMonitoringTransition> trgroup;
  TransitionGroupOpenMS<MSSpectrum, ReactionMonitoringTransition> default_transition_group_openms(trgroup);

} //end namespace OpenMS

