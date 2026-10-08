// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConversionHelper.h>

namespace OpenMS
{
  namespace
  {
    /// The native counterpart of the legacy conversion: the identifications of the converted features are kept and
    /// marked with the map index, those of the other features and of subordinates are dropped, unassigned ones are kept.
    void convertIdentifications(UInt64 map_index, const FeatureMap& input, Size n, ConsensusMap& output)
    {
      using ID = IdentificationData;
      auto data = input.getIdentificationData();
      if (data.empty()) return;
      std::set<ID::QueryReference> linked, kept;
      const auto collect = [&](const auto& self, const Feature& feature) -> void {
        const auto queries = feature.getLinkedIDQueries(data);
        linked.insert(queries.begin(), queries.end());
        for (const auto& subordinate : feature.getSubordinates())
          self(self, subordinate);
      };
      for (Size i = 0; i < input.size(); ++i)
      {
        collect(collect, input[i]);
        if (i < n)
        {
          const auto queries = input[i].getLinkedIDQueries(data);
          kept.insert(queries.begin(), queries.end());
        }
      }
      for (const auto& current : data.getRuns())
      {
        auto& run = data.getRun(current.getIdentifier());
        run.eraseIdentifications([&](const ID::Identification& query) {
          const ID::QueryReference reference {run.getUuid(), query.getId()};
          return linked.contains(reference) && ! kept.contains(reference);
        });
        for (const auto& source : run.getSources())
          for (const auto& query : source.identifications)
          {
            if (! kept.contains({run.getUuid(), query.getId()})) continue;
            ID::Observation observation = query;
            observation.setMetaValue("map_index", map_index);
            run.replaceObservation(query.getId(), observation);
          }
      }
      output.getIdentificationData() = std::move(data);
    }
  } // namespace

  void MapConversion::convert(UInt64 const input_map_index,
                              PeakMap& input_map,
                              ConsensusMap& output_map,
                              Size n)
  {
    output_map.clear(true);

    // see @todo above
    output_map.setUniqueId();

    input_map.updateRanges();
    std::vector<Peak2D> tmp;
    tmp.reserve(input_map.getSize()); // an upper bound only, see below

    // TODO Avoid tripling the memory consumption by this call
    input_map.get2DData(tmp);

    // Clamp n only now, against the number of peaks actually collected.
    // input_map.getSize() is the wrong bound: it counts the peaks of every
    // spectrum at every MS level plus all chromatogram points, whereas
    // get2DData() collects MS1 peaks only. With n > tmp.size() the middle
    // iterator of the partial_sort and the copy loop below would run past
    // the end of tmp (out-of-bounds reads and writes, consensus features
    // built from garbage).
    if (n > tmp.size())
    {
      n = tmp.size();
    }
    output_map.reserve(n);

    // most intense first; equal intensities are ordered by RT, then m/z, so that the selected peaks and their
    // order do not depend on the standard library's partial_sort implementation
    std::partial_sort(tmp.begin(),
                      tmp.begin() + n,
                      tmp.end(),
                      [](const Peak2D& left, const Peak2D& right)
                      {
                        if (left.getIntensity() != right.getIntensity()) return left.getIntensity() > right.getIntensity();
                        if (left.getRT() != right.getRT()) return left.getRT() < right.getRT();
                        return left.getMZ() < right.getMZ();
                      });

    for (Size element_index = 0; element_index < n; ++element_index)
    {
      output_map.push_back(ConsensusFeature(input_map_index,
                                            tmp[element_index],
                                            element_index));
    }

    output_map.getColumnHeaders()[input_map_index].size = n;
    output_map.updateRanges();
  }

  void MapConversion::convert(ConsensusMap const& input_map,
                              const bool keep_uids,
                              FeatureMap& output_map)
  {
    output_map.clear(true);
    output_map.resize(input_map.size());
    output_map.DocumentIdentifier::operator=(input_map);

    if (keep_uids)
    {
      output_map.UniqueIdInterface::operator=(input_map);
    }
    else
    {
      output_map.setUniqueId();
    }
    output_map.setProteinIdentifications(input_map.getProteinIdentifications());
    output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());
    output_map.getIdentificationData() = input_map.getIdentificationData();

    for (Size i = 0; i < input_map.size(); ++i)
    {
      Feature& f = output_map[i];
      const ConsensusFeature& c = input_map[i];
      f.BaseFeature::operator=(c);
      if (!keep_uids)
      {
        f.setUniqueId();
      }
    }

    output_map.updateRanges();
  }

  void MapConversion::convert(UInt64 const input_map_index,
                              FeatureMap const& input_map,
                              ConsensusMap& output_map,
                              Size n)
  {
    if (n > input_map.size())
    {
      n = input_map.size();
    }

    output_map.clear(true);
    output_map.reserve(n);

    // An arguable design decision, see above.
    output_map.setUniqueId(input_map.getUniqueId());

    for (UInt64 element_index = 0; element_index < n; ++element_index)
    {
      output_map.push_back(ConsensusFeature(input_map_index, input_map[element_index]));
    }
    output_map.getColumnHeaders()[input_map_index].size = static_cast<Size>(input_map.size());
    output_map.setProteinIdentifications(input_map.getProteinIdentifications());
    output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());
    convertIdentifications(input_map_index, input_map, n, output_map);
    output_map.updateRanges();
  }

} // namespace OpenMS
