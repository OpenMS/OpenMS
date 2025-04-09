// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>
#include <OpenMS/MATH/MathFunctions.h>

using namespace std;

namespace OpenMS
{
  // return map of ms2 to feature and a vector of unassigned ms2
  FeatureMapping::FeatureToMs2Indices FeatureMapping::assignMS2IndexToFeature(const OpenMS::MSExperiment& spectra,
                                                                              const FeatureMappingInfo& fm_info,
                                                                              const double& precursor_mz_tolerance,
                                                                              const double& precursor_rt_tolerance,
                                                                              bool ppm)
  {
    std::map<const BaseFeature*, std::vector<size_t>>  assigned_ms2;
    vector<size_t> unassigned_ms2;

    // map precursors to closest feature and retrieve annotated metadata (if possible)
    for (size_t index = 0; index != spectra.size(); ++index)
    {
      if (spectra[index].getMSLevel() != 2) { continue; }

      // get precursor meta data (m/z, rt)
      const vector<Precursor> & pcs = spectra[index].getPrecursors();

      if (!pcs.empty())
      {
        const double mz = pcs[0].getMZ();
        const double rt = spectra[index].getRT();

        // query features in tolerance window
        vector<Size> matches;

        // get mz tolerance window
        std::pair<double,double> mz_tolerance_window = Math::getTolWindow(mz, precursor_mz_tolerance, ppm);
        fm_info.kd_tree.queryRegion(rt - precursor_rt_tolerance, rt + precursor_rt_tolerance, mz_tolerance_window.first, mz_tolerance_window.second, matches, true);

        // no precursor matches the feature information found
        if (matches.empty())
        {
          unassigned_ms2.push_back(index);
          continue;
        }

        // in the case of multiple features in tolerance window, select the one closest in m/z to the precursor
        Size min_distance_feature_index(0);
        double min_distance(1e11);
        for (auto const & k_idx : matches)
        {
          const double f_mz = fm_info.kd_tree.mz(k_idx);
          const double distance = fabs(f_mz - mz);
          if (distance < min_distance)
          {
            min_distance = distance;
            min_distance_feature_index = k_idx;
          }
        }
        const BaseFeature* min_distance_feature = fm_info.kd_tree.feature(min_distance_feature_index);
        assigned_ms2[min_distance_feature].push_back(index);
      }
    }
    FeatureMapping::FeatureToMs2Indices feature_mapping;
    feature_mapping.assignedMS2 = assigned_ms2;
    feature_mapping.unassignedMS2 = unassigned_ms2;
    return feature_mapping;
  }
} // namespace OpenMS
