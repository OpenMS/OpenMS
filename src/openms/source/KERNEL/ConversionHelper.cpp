// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <numeric>

namespace OpenMS
{
  namespace
  {
    /// The identifications of the converted features are kept and
    /// marked with the map index, those of the other features and of subordinates are dropped, unassigned ones are kept.
    void convertIdentifications(UInt64 map_index, const FeatureMap& input, const std::vector<Size>& selected, ConsensusMap& output)
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
      for (const Feature& feature : input)
      {
        collect(collect, feature);
      }
      for (Size i : selected)
      {
        const auto queries = input[i].getLinkedIDQueries(data);
        kept.insert(queries.begin(), queries.end());
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
    std::vector<Size> selected(n);
    std::iota(selected.begin(), selected.end(), Size(0));
    convert(input_map_index, input_map, output_map, selected);
  }

  void MapConversion::convert(UInt64 const input_map_index,
                              FeatureMap const& input_map,
                              ConsensusMap& output_map,
                              const std::vector<Size>& selected)
  {
    for (Size index : selected)
    {
      if (index >= input_map.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, input_map.size());
      }
    }

    output_map.clear(true);
    output_map.reserve(selected.size());

    // An arguable design decision, see above.
    output_map.setUniqueId(input_map.getUniqueId());

    for (Size index : selected)
    {
      output_map.push_back(ConsensusFeature(input_map_index, input_map[index]));
    }
    output_map.getColumnHeaders()[input_map_index].size = static_cast<Size>(input_map.size());
    convertIdentifications(input_map_index, input_map, selected, output_map);
    output_map.updateRanges();
  }

} // namespace OpenMS
