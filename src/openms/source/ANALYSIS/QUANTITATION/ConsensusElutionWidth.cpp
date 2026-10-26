// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/ConsensusElutionWidth.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <cmath>
#include <unordered_map>

namespace OpenMS
{
  namespace
  {
    const double SQRT_2PI = std::sqrt(2.0 * Constants::PI);

    /// height and sigma of a valid Gaussian fit, else false
    bool gaussFit(const Feature& f, double& height, double& sigma, double& center)
    {
      if (!f.metaValueExists("model_height") || !f.metaValueExists("model_Gauss_sigma") || !f.metaValueExists("model_center")) return false;
      if (f.metaValueExists("model_status"))
      {
        const std::string status = f.getMetaValue("model_status").toString();
        if (status.empty() || status[0] != '0') return false;
      }
      height = f.getMetaValue("model_height");
      sigma = f.getMetaValue("model_Gauss_sigma");
      center = f.getMetaValue("model_center");
      return std::isfinite(height) && height > 0.0 && std::isfinite(sigma) && sigma > 0.0 && std::isfinite(center);
    }
  }

  double ConsensusElutionWidth::localSlope(const TransformationDescription& trafo, double rt)
  {
    if (trafo.getModelType() == "none" || trafo.getModelType() == "identity") return 1.0;
    const double h = 1.0;
    const double slope = (trafo.apply(rt + h) - trafo.apply(rt - h)) / (2.0 * h);
    return (std::isfinite(slope) && slope > 0.0) ? slope : 1.0;
  }

  ConsensusElutionWidth::Stats ConsensusElutionWidth::apply(ConsensusMap& consensus,
                                                            const std::vector<FeatureMap>& feature_maps,
                                                            const std::vector<TransformationDescription>& transformations)
  {
    Stats stats;
    // unique id -> feature, per map
    std::vector<std::unordered_map<UInt64, const Feature*>> by_id(feature_maps.size());
    for (Size m = 0; m < feature_maps.size(); ++m)
    {
      by_id[m].reserve(feature_maps[m].size());
      for (const Feature& f : feature_maps[m]) by_id[m].emplace(f.getUniqueId(), &f);
    }
    auto slopeOf = [&](Size map_index, double rt)
    {
      return map_index < transformations.size() ? localSlope(transformations[map_index], rt) : 1.0;
    };

    std::vector<double> widths;
    for (ConsensusFeature& cf : consensus)
    {
      struct Fit { double height, sigma, center, slope; bool ok; };
      std::vector<Fit> fits;
      fits.reserve(cf.size());
      widths.clear();
      for (const FeatureHandle& fh : cf.getFeatures())
      {
        Fit fit{0.0, 0.0, 0.0, 1.0, false};
        const Size m = fh.getMapIndex();
        if (m < by_id.size())
        {
          auto it = by_id[m].find(fh.getUniqueId());
          if (it != by_id[m].end() && gaussFit(*it->second, fit.height, fit.sigma, fit.center))
          {
            fit.ok = true;
            fit.slope = slopeOf(m, fit.center);
            if (!it->second->getPeptideIdentifications().empty()) widths.push_back(fit.sigma * fit.slope);
          }
        }
        fits.push_back(fit);
      }
      if (widths.empty())
      {
        stats.handles_kept += cf.size();
        continue;
      }
      const double w_c = Math::median(widths.begin(), widths.end());
      ++stats.consensus_features;

      // ConsensusFeature::getFeatures() is a set ordered by map index and unique id; rebuild it with the new intensities
      ConsensusFeature::HandleSetType handles;
      Size i = 0;
      double total = 0.0;
      for (const FeatureHandle& fh : cf.getFeatures())
      {
        FeatureHandle copy = fh;
        const Fit& fit = fits[i++];
        if (fit.ok)
        {
          copy.setIntensity(static_cast<float>(fit.height * (w_c / fit.slope) * SQRT_2PI));
          ++stats.handles_updated;
        }
        else
        {
          ++stats.handles_kept;
        }
        total += copy.getIntensity();
        handles.insert(copy);
      }
      cf.setFeatures(std::move(handles));
      // as ConsensusFeature::computeConsensus(): the mean of the handles' intensities
      cf.setIntensity(static_cast<float>(total / static_cast<double>(cf.size())));
    }
    return stats;
  }
}
