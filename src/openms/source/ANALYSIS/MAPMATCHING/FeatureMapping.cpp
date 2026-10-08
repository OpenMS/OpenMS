// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

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
    for (const auto& map : fm_info.feature_maps)
      if (IdentificationDataConverter::hasPeptideIdentifications(map))
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Feature maps with peptide identifications: move them into identification data first");
    std::map<const BaseFeature*, std::vector<size_t>>  assigned_ms2;
    std::map<const BaseFeature*, const IdentificationData*> identification_data;
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
        const Size map_index = fm_info.kd_tree.mapIndex(min_distance_feature_index);
        if (map_index < fm_info.feature_maps.size())
          identification_data[min_distance_feature] = &fm_info.feature_maps[map_index].getIdentificationData();
      }
    }
    FeatureMapping::FeatureToMs2Indices feature_mapping;
    feature_mapping.assignedMS2 = assigned_ms2;
    feature_mapping.unassignedMS2 = unassigned_ms2;
    feature_mapping.identification_data = std::move(identification_data);
    return feature_mapping;
  }

  std::vector<const IdentificationData::Match*> FeatureMapping::FeatureToMs2Indices::getFirstLinkedMatches(const BaseFeature* feature) const
  {
    const auto data = identification_data.find(feature);
    if (data == identification_data.end()) return {};
    auto linked = feature->getLinkedIdentifications(*data->second);
    if (linked.empty()) return {};
    return std::move(linked.front().matches);
  }
} // namespace OpenMS
