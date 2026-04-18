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
    features_ = std::make_shared<FeatureWrapperList_>();
    features_->reserve(ids.size());
    feature_index_.reserve(ids.size());
    for (std::vector<String>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      feature_index_.emplace_back(*it, features_->size());
      features_->emplace_back(mrmfeature.getFeature(*it));
    }

    std::vector<String> p_ids;
    mrmfeature.getPrecursorFeatureIDs(p_ids);
    precursor_features_ = std::make_shared<FeatureWrapperList_>();
    precursor_features_->reserve(p_ids.size());
    precursor_feature_index_.reserve(p_ids.size());
    for (std::vector<String>::iterator it = p_ids.begin(); it != p_ids.end(); ++it)
    {
      precursor_feature_index_.emplace_back(*it, precursor_features_->size());
      precursor_features_->emplace_back(mrmfeature.getPrecursorFeature(*it));
    }
  }

  MRMFeatureOpenMS::MRMFeatureOpenMS(MRMFeature& mrmfeature,
                                     const std::vector<std::string>& feature_ids,
                                     const std::vector<std::string>& precursor_feature_ids) :
    mrmfeature_(mrmfeature)
  {
    features_ = std::make_shared<FeatureWrapperList_>();
    features_->reserve(feature_ids.size());
    feature_index_.reserve(feature_ids.size());
    for (const std::string& id : feature_ids)
    {
      feature_index_.emplace_back(id, features_->size());
      features_->emplace_back(mrmfeature.getFeature(String(id)));
    }

    precursor_features_ = std::make_shared<FeatureWrapperList_>();
    precursor_features_->reserve(precursor_feature_ids.size());
    precursor_feature_index_.reserve(precursor_feature_ids.size());
    for (const std::string& id : precursor_feature_ids)
    {
      precursor_feature_index_.emplace_back(id, precursor_features_->size());
      precursor_features_->emplace_back(mrmfeature.getPrecursorFeature(String(id)));
    }
  }

  MRMFeatureOpenMS::MRMFeatureOpenMS(MRMFeature& mrmfeature,
                                     const std::vector<std::string>& feature_ids,
                                     const std::vector<std::string>& precursor_feature_ids,
                                     const std::vector<String>& feature_lookup_ids,
                                     const std::vector<String>& precursor_feature_lookup_ids) :
    mrmfeature_(mrmfeature),
    feature_ids_(&feature_ids),
    precursor_feature_ids_(&precursor_feature_ids)
  {
    OPENMS_PRECONDITION(feature_ids.size() == feature_lookup_ids.size(),
                        "Feature id cache needs to match feature lookup id cache.");
    OPENMS_PRECONDITION(precursor_feature_ids.size() == precursor_feature_lookup_ids.size(),
                        "Precursor feature id cache needs to match precursor feature lookup id cache.");

    features_ = std::make_shared<FeatureWrapperList_>();
    features_->reserve(feature_lookup_ids.size());
    for (const String& id : feature_lookup_ids)
    {
      features_->emplace_back(mrmfeature.getFeature(id));
    }

    precursor_features_ = std::make_shared<FeatureWrapperList_>();
    precursor_features_->reserve(precursor_feature_lookup_ids.size());
    for (const String& id : precursor_feature_lookup_ids)
    {
      precursor_features_->emplace_back(mrmfeature.getPrecursorFeature(id));
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
    const ConvexHull2D::PointArrayType& data_points = feature_->getConvexHulls()[0].getHullPoints();
    rt.reserve(rt.size() + data_points.size());
    for (ConvexHull2D::PointArrayType::const_iterator it = data_points.begin(); it != data_points.end(); ++it)
    {
      rt.push_back(it->getX());
    }
  }

  void FeatureOpenMS::getIntensity(std::vector<double>& intens) const
  {
    OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "There needs to exactly one convex hull per feature.");
    const ConvexHull2D::PointArrayType& data_points = feature_->getConvexHulls()[0].getHullPoints();
    intens.reserve(intens.size() + data_points.size());
    for (ConvexHull2D::PointArrayType::const_iterator it = data_points.begin(); it != data_points.end(); ++it)
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

  std::shared_ptr<OpenSwath::IFeature> MRMFeatureOpenMS::getFeature(std::string nativeID)
  {
    const std::size_t index = findFeatureIndex_(nativeID);
    return std::shared_ptr<OpenSwath::IFeature>(features_, &(*features_)[index]);
  }

  std::shared_ptr<OpenSwath::IFeature> MRMFeatureOpenMS::getPrecursorFeature(std::string nativeID)
  {
    const std::size_t index = findPrecursorFeatureIndex_(nativeID);
    return std::shared_ptr<OpenSwath::IFeature>(precursor_features_, &(*precursor_features_)[index]);
  }

  std::vector<std::string> MRMFeatureOpenMS::getNativeIDs() const
  {
    if (feature_ids_ != nullptr)
    {
      return *feature_ids_;
    }
    std::vector<std::string> v;
    v.reserve(feature_index_.size());
    for (const auto& entry : feature_index_)
    {
      v.push_back(entry.first);
    }
    return v;
  }

  std::vector<std::string> MRMFeatureOpenMS::getPrecursorIDs() const
  {
    if (precursor_feature_ids_ != nullptr)
    {
      return *precursor_feature_ids_;
    }
    std::vector<std::string> v;
    v.reserve(precursor_feature_index_.size());
    for (const auto& entry : precursor_feature_index_)
    {
      v.push_back(entry.first);
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

  double MRMFeatureOpenMS::getMetaValue(std::string name) const
  {
    return mrmfeature_.getMetaValue(name);
  }

  size_t MRMFeatureOpenMS::size() const
  {
    return feature_ids_ != nullptr ? feature_ids_->size() : feature_index_.size();
  }

  void MRMFeatureOpenMS::getFeatureIntensities(const std::vector<std::string>& native_ids, std::vector<std::vector<double>>& intensities) const
  {
    intensities.resize(native_ids.size());
    if (feature_ids_ == &native_ids)
    {
      for (std::size_t i = 0; i < intensities.size(); ++i)
      {
        intensities[i].clear();
        (*features_)[i].getIntensity(intensities[i]);
      }
      return;
    }
    for (std::size_t i = 0; i < intensities.size(); ++i)
    {
      intensities[i].clear();
      (*features_)[findFeatureIndex_(native_ids[i])].getIntensity(intensities[i]);
    }
  }

  void MRMFeatureOpenMS::getPrecursorFeatureIntensities(const std::vector<std::string>& native_ids, std::vector<std::vector<double>>& intensities) const
  {
    intensities.resize(native_ids.size());
    if (precursor_feature_ids_ == &native_ids)
    {
      for (std::size_t i = 0; i < intensities.size(); ++i)
      {
        intensities[i].clear();
        (*precursor_features_)[i].getIntensity(intensities[i]);
      }
      return;
    }
    for (std::size_t i = 0; i < intensities.size(); ++i)
    {
      intensities[i].clear();
      (*precursor_features_)[findPrecursorFeatureIndex_(native_ids[i])].getIntensity(intensities[i]);
    }
  }

  float MRMFeatureOpenMS::getFeatureIntensity(const std::string& native_id) const
  {
    return (*features_)[findFeatureIndex_(native_id)].getIntensity();
  }

  float MRMFeatureOpenMS::getFeatureIntensity(const std::string& native_id, std::size_t expected_index) const
  {
    if (feature_ids_ != nullptr && expected_index < feature_ids_->size() && (*feature_ids_)[expected_index] == native_id)
    {
      return (*features_)[expected_index].getIntensity();
    }
    return getFeatureIntensity(native_id);
  }

  std::size_t MRMFeatureOpenMS::findFeatureIndex_(const std::string& native_id) const
  {
    if (feature_ids_ != nullptr)
    {
      for (std::size_t i = 0; i < feature_ids_->size(); ++i)
      {
        if ((*feature_ids_)[i] == native_id)
        {
          return i;
        }
      }
    }
    else
    {
      for (const auto& entry : feature_index_)
      {
        if (entry.first == native_id)
        {
          return entry.second;
        }
      }
    }
    OPENMS_PRECONDITION(false, "Feature needs to exist");
    return 0;
  }

  std::size_t MRMFeatureOpenMS::findPrecursorFeatureIndex_(const std::string& native_id) const
  {
    if (precursor_feature_ids_ != nullptr)
    {
      for (std::size_t i = 0; i < precursor_feature_ids_->size(); ++i)
      {
        if ((*precursor_feature_ids_)[i] == native_id)
        {
          return i;
        }
      }
    }
    else
    {
      for (const auto& entry : precursor_feature_index_)
      {
        if (entry.first == native_id)
        {
          return entry.second;
        }
      }
    }
    OPENMS_PRECONDITION(false, "Precursor feature needs to exist");
    return 0;
  }

  // default instances
  MSSpectrum spec;
  MSChromatogram chrom;
  SignalToNoiseOpenMS< MSSpectrum> spec_signal_to_noise_openms(spec, 1.0, 3, true);
  SignalToNoiseOpenMS< MSChromatogram > chrom_signal_to_noise_openms(chrom, 1.0, 3, true);

  MRMTransitionGroup<MSSpectrum, ReactionMonitoringTransition> trgroup;
  TransitionGroupOpenMS<MSSpectrum, ReactionMonitoringTransition> default_transition_group_openms(trgroup);

} //end namespace OpenMS
