// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

namespace OpenMS
{
  namespace
  {
    void appendRT_(const Feature& feature, std::vector<double>& rt)
    {
      const std::vector<ConvexHull2D>& convex_hulls = feature.getConvexHulls();
      if (convex_hulls.size() != 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "There needs to be exactly one convex hull per feature.");
      }
      const ConvexHull2D::PointArrayType& data_points = convex_hulls[0].getHullPoints();
      rt.reserve(rt.size() + data_points.size());
      for (ConvexHull2D::PointArrayType::const_iterator it = data_points.begin(); it != data_points.end(); ++it)
      {
        rt.push_back(it->getX());
      }
    }

    void appendIntensity_(const Feature& feature, std::vector<double>& intens)
    {
      const std::vector<ConvexHull2D>& convex_hulls = feature.getConvexHulls();
      if (convex_hulls.size() != 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "There needs to be exactly one convex hull per feature.");
      }
      const ConvexHull2D::PointArrayType& data_points = convex_hulls[0].getHullPoints();
      intens.reserve(intens.size() + data_points.size());
      for (ConvexHull2D::PointArrayType::const_iterator it = data_points.begin(); it != data_points.end(); ++it)
      {
        intens.push_back(it->getY());
      }
    }
  }

  struct MRMFeatureOpenMS::FeatureLane_
  {
    using FeatureWrapperList_ = std::vector<FeatureOpenMS>;
    using FeaturePointerList_ = std::vector<const Feature*>;
    using FeatureIndexList_ = std::vector<std::pair<std::string, Size>>;

    std::shared_ptr<FeatureWrapperList_> wrappers;
    const std::vector<Feature>* direct_features = nullptr;
    FeaturePointerList_ feature_ptrs;
    FeatureIndexList_ feature_index;
    std::vector<std::string> ids_storage;
    const std::vector<std::string>* ids = nullptr;

    template <typename Getter>
    void initializeFromLookupIds(const std::vector<String>& lookup_ids, Getter&& getter)
    {
      feature_ptrs.reserve(lookup_ids.size());
      feature_index.reserve(lookup_ids.size());
      for (const String& id : lookup_ids)
      {
        feature_index.emplace_back(id, feature_ptrs.size());
        feature_ptrs.push_back(&getter(id));
      }
    }

    template <typename Getter>
    void initializeFromStringIds(const std::vector<std::string>& lookup_ids, Getter&& getter)
    {
      feature_ptrs.reserve(lookup_ids.size());
      feature_index.reserve(lookup_ids.size());
      for (const std::string& id : lookup_ids)
      {
        feature_index.emplace_back(id, feature_ptrs.size());
        feature_ptrs.push_back(&getter(String(id)));
      }
    }

    void initializeAlignedDirect(const std::vector<Feature>& features, const std::vector<std::string>& aligned_ids)
    {
      direct_features = &features;
      ids_storage = aligned_ids;
      ids = &ids_storage;
    }

    void setAlignedIds(const std::vector<std::string>& aligned_ids)
    {
      ids_storage = aligned_ids;
      ids = &ids_storage;
    }

    Size size() const
    {
      return ids != nullptr ? ids->size() : feature_ptrs.size();
    }

    std::vector<std::string> nativeIds() const
    {
      if (ids != nullptr)
      {
        return *ids;
      }
      std::vector<std::string> result;
      result.reserve(feature_index.size());
      for (const auto& entry : feature_index)
      {
        result.push_back(entry.first);
      }
      return result;
    }

    Size findIndex(const std::string& native_id) const
    {
      if (ids != nullptr)
      {
        for (Size i = 0; i < ids->size(); ++i)
        {
          if ((*ids)[i] == native_id)
          {
            return i;
          }
        }
      }
      else
      {
        for (const auto& entry : feature_index)
        {
          if (entry.first == native_id)
          {
            return entry.second;
          }
        }
      }
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, native_id);
    }

    const Feature& featureByIndex(Size index) const
    {
      return direct_features != nullptr ? (*direct_features)[index] : *feature_ptrs[index];
    }

    std::shared_ptr<OpenSwath::IFeature> wrapperByIndex(Size index)
    {
      ensureWrappers_();
      return std::shared_ptr<OpenSwath::IFeature>(wrappers, &(*wrappers)[index]);
    }

    void appendIntensitiesByIds(const std::vector<std::string>& native_ids, std::vector<std::vector<double>>& intensities) const
    {
      intensities.resize(native_ids.size());
      if (ids != nullptr && native_ids == *ids)
      {
        for (Size i = 0; i < intensities.size(); ++i)
        {
          intensities[i].clear();
          appendIntensity_(featureByIndex(i), intensities[i]);
        }
        return;
      }
      for (Size i = 0; i < intensities.size(); ++i)
      {
        intensities[i].clear();
        appendIntensity_(featureByIndex(findIndex(native_ids[i])), intensities[i]);
      }
    }

    float intensity(const std::string& native_id, Size expected_index) const
    {
      if (ids != nullptr && expected_index < ids->size() && (*ids)[expected_index] == native_id)
      {
        return featureByIndex(expected_index).getIntensity();
      }
      return featureByIndex(findIndex(native_id)).getIntensity();
    }

  private:
    void ensureWrappers_()
    {
      if (wrappers != nullptr)
      {
        return;
      }
      auto tmp = std::make_shared<FeatureWrapperList_>();
      if (direct_features != nullptr)
      {
        tmp->reserve(direct_features->size());
        for (const Feature& feature : *direct_features)
        {
          tmp->emplace_back(feature);
        }
      }
      else
      {
        tmp->reserve(feature_ptrs.size());
        for (const Feature* feature : feature_ptrs)
        {
          tmp->emplace_back(*feature);
        }
      }
      wrappers = std::move(tmp);
    }
  };

  MRMFeatureOpenMS::MRMFeatureOpenMS(MRMFeature& mrmfeature) :
    mrmfeature_(mrmfeature),
    feature_lane_(std::make_unique<FeatureLane_>()),
    precursor_feature_lane_(std::make_unique<FeatureLane_>())
  {
    std::vector<String> ids;
    mrmfeature.getFeatureIDs(ids);
    feature_lane_->initializeFromLookupIds(ids, [&](const String& id) -> const Feature& { return mrmfeature.getFeature(id); });

    std::vector<String> p_ids;
    mrmfeature.getPrecursorFeatureIDs(p_ids);
    precursor_feature_lane_->initializeFromLookupIds(p_ids, [&](const String& id) -> const Feature& { return mrmfeature.getPrecursorFeature(id); });
  }

  MRMFeatureOpenMS::MRMFeatureOpenMS(MRMFeature& mrmfeature,
                                     const std::vector<std::string>& feature_ids,
                                     const std::vector<std::string>& precursor_feature_ids) :
    mrmfeature_(mrmfeature),
    feature_lane_(std::make_unique<FeatureLane_>()),
    precursor_feature_lane_(std::make_unique<FeatureLane_>())
  {
    feature_lane_->initializeFromStringIds(feature_ids, [&](const String& id) -> const Feature& { return mrmfeature.getFeature(id); });
    precursor_feature_lane_->initializeFromStringIds(precursor_feature_ids, [&](const String& id) -> const Feature& { return mrmfeature.getPrecursorFeature(id); });
  }

  MRMFeatureOpenMS::MRMFeatureOpenMS(MRMFeature& mrmfeature,
                                     const std::vector<std::string>& feature_ids,
                                     const std::vector<std::string>& precursor_feature_ids,
                                     const std::vector<String>& feature_lookup_ids,
                                     const std::vector<String>& precursor_feature_lookup_ids) :
    mrmfeature_(mrmfeature),
    feature_lane_(std::make_unique<FeatureLane_>()),
    precursor_feature_lane_(std::make_unique<FeatureLane_>())
  {
    if (feature_ids.size() != feature_lookup_ids.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Feature id cache needs to match feature lookup id cache.");
    }
    if (precursor_feature_ids.size() != precursor_feature_lookup_ids.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Precursor feature id cache needs to match precursor feature lookup id cache.");
    }

    const std::vector<Feature>& fragment_features = mrmfeature.getFeatures();
    if (fragment_features.size() != feature_lookup_ids.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Feature storage needs to match feature lookup id cache.");
    }
    feature_lane_->initializeAlignedDirect(fragment_features, feature_ids);
    precursor_feature_lane_->initializeFromLookupIds(precursor_feature_lookup_ids, [&](const String& id) -> const Feature& { return mrmfeature.getPrecursorFeature(id); });
    precursor_feature_lane_->setAlignedIds(precursor_feature_ids);
  }

  FeatureOpenMS::FeatureOpenMS(const Feature& feature) :
    feature_(&feature)
  {
  }

  FeatureOpenMS::~FeatureOpenMS() = default;

  void FeatureOpenMS::getRT(std::vector<double>& rt) const
  {
    appendRT_(*feature_, rt);
  }

  void FeatureOpenMS::getIntensity(std::vector<double>& intens) const
  {
    appendIntensity_(*feature_, intens);
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
    return feature_lane_->wrapperByIndex(feature_lane_->findIndex(nativeID));
  }

  std::shared_ptr<OpenSwath::IFeature> MRMFeatureOpenMS::getPrecursorFeature(std::string nativeID)
  {
    return precursor_feature_lane_->wrapperByIndex(precursor_feature_lane_->findIndex(nativeID));
  }

  std::vector<std::string> MRMFeatureOpenMS::getNativeIDs() const
  {
    return feature_lane_->nativeIds();
  }

  std::vector<std::string> MRMFeatureOpenMS::getPrecursorIDs() const
  {
    return precursor_feature_lane_->nativeIds();
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
    return feature_lane_->size();
  }

  void MRMFeatureOpenMS::getFeatureIntensities(const std::vector<std::string>& native_ids, std::vector<std::vector<double>>& intensities) const
  {
    feature_lane_->appendIntensitiesByIds(native_ids, intensities);
  }

  void MRMFeatureOpenMS::getPrecursorFeatureIntensities(const std::vector<std::string>& native_ids, std::vector<std::vector<double>>& intensities) const
  {
    precursor_feature_lane_->appendIntensitiesByIds(native_ids, intensities);
  }

  float MRMFeatureOpenMS::getFeatureIntensity(const std::string& native_id, Size expected_index) const
  {
    return feature_lane_->intensity(native_id, expected_index);
  }

  // default instances
  MSSpectrum spec;
  MSChromatogram chrom;
  SignalToNoiseOpenMS< MSSpectrum> spec_signal_to_noise_openms(spec, 1.0, 3, true);
  SignalToNoiseOpenMS< MSChromatogram > chrom_signal_to_noise_openms(chrom, 1.0, 3, true);

  MRMTransitionGroup<MSSpectrum, ReactionMonitoringTransition> trgroup;
  TransitionGroupOpenMS<MSSpectrum, ReactionMonitoringTransition> default_transition_group_openms(trgroup);

} //end namespace OpenMS
