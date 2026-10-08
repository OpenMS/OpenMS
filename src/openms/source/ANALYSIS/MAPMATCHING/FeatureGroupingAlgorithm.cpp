// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConversionHelper.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithm::FeatureGroupingAlgorithm() :
    DefaultParamHandler("FeatureGroupingAlgorithm")
  {
  }

  void FeatureGroupingAlgorithm::group(const vector<ConsensusMap>& maps, ConsensusMap& out)
  {
    OPENMS_LOG_WARN << "FeatureGroupingAlgorithm::group() does not support ConsensusMaps directly. Converting to FeatureMaps." << endl;

    vector<FeatureMap> maps_f;
    for (Size i = 0; i < maps.size(); ++i)
    {
      FeatureMap fm;
      MapConversion::convert(maps[i], true, fm);
      maps_f.push_back(fm);
    }
    // call FeatureMap version of group()
    group(maps_f, out);
  }

  void FeatureGroupingAlgorithm::transferSubelements(const vector<ConsensusMap>& maps, ConsensusMap& out) const
  {
    // accumulate file descriptions from the input maps:
    // cout << "Updating file descriptions..." << endl;
    out.getColumnHeaders().clear();
    // mapping: (input file index / map index assigned by the linkers, old map index) -> new map index
    map<pair<Size, UInt64>, Size> mapid_table;
    for (Size i = 0; i < maps.size(); ++i)
    {
      const ConsensusMap& consensus = maps[i];
      for (ConsensusMap::ColumnHeaders::const_iterator desc_it = consensus.getColumnHeaders().begin(); desc_it != consensus.getColumnHeaders().end(); ++desc_it)
      {
        Size counter = mapid_table.size();
        mapid_table[make_pair(i, desc_it->first)] = counter;
        out.getColumnHeaders()[counter] = desc_it->second;
      }
    }

    // look-up table: input map -> unique ID -> consensus feature
    // cout << "Creating look-up table..." << endl;
    vector<map<UInt64, ConsensusMap::ConstIterator> > feat_lookup(maps.size());
    for (Size i = 0; i < maps.size(); ++i)
    {
      const ConsensusMap& consensus = maps[i];
      for (ConsensusMap::ConstIterator feat_it = consensus.begin();
           feat_it != consensus.end(); ++feat_it)
      {
        // do NOT use "id_lookup[i][feat_it->getUniqueId()] = feat_it;" here as
        // you will get "attempt to copy-construct an iterator from a singular
        // iterator" in STL debug mode:
        feat_lookup[i].insert(make_pair(feat_it->getUniqueId(), feat_it));
      }
    }
    // adjust the consensus features:
    // cout << "Adjusting consensus features..." << endl;
    for (ConsensusMap::iterator cons_it = out.begin(); cons_it != out.end(); ++cons_it)
    {
      ConsensusFeature adjusted = ConsensusFeature(
        static_cast<BaseFeature>(*cons_it)); // remove sub-features
      for (ConsensusFeature::HandleSetType::const_iterator sub_it = cons_it->getFeatures().begin(); sub_it != cons_it->getFeatures().end(); ++sub_it)
      {
        UInt64 id = sub_it->getUniqueId();
        Size map_index = sub_it->getMapIndex();
        ConsensusMap::ConstIterator origin = feat_lookup[map_index][id];
        for (ConsensusFeature::HandleSetType::const_iterator handle_it = origin->getFeatures().begin(); handle_it != origin->getFeatures().end(); ++handle_it)
        {
          FeatureHandle handle = *handle_it;
          Size new_id = mapid_table[make_pair(map_index, handle.getMapIndex())];
          handle.setMapIndex(new_id);
          adjusted.insert(handle);
        }
      }
      *cons_it = adjusted;

      for (auto& id : cons_it->getPeptideIdentifications())
      {
        // if old_map_index is not present, there was no map_index in the beginning,
        // therefore the newly assigned map_index cannot be "corrected"
        // -> remove the MetaValue to be consistent.
        if (id.metaValueExists("old_map_index"))
        {
          Size old_map_index = (Size)id.getMetaValue("old_map_index");
          Size file_index = (Size)id.getMetaValue("map_index");
          Size new_idx = mapid_table[make_pair(file_index, old_map_index)];
          id.setMetaValue("map_index", new_idx);
          id.removeMetaValue("old_map_index");
        }
        else
        {
          id.removeMetaValue("map_index");
        }
      }
    }
    for (auto& id : out.getUnassignedPeptideIdentifications())
    {
      // if old_map_index is not present, there was no map_index in the beginning,
      // therefore the newly assigned map_index cannot be "corrected"
      // -> remove the MetaValue to be consistent.
      if (id.metaValueExists("old_map_index"))
      {
        Size old_map_index = (Size)id.getMetaValue("old_map_index");
        Size file_index = (Size)id.getMetaValue("map_index");
        Size new_idx = mapid_table[make_pair(file_index, old_map_index)];
        id.setMetaValue("map_index", new_idx);
        id.removeMetaValue("old_map_index");
      }
      else
      {
        id.removeMetaValue("map_index");
      }
    }
    // the same for the identifications in identification data:
    auto& data = out.getIdentificationData();
    for (const auto& current : data.getRuns())
    {
      auto& run = data.getRun(current.getIdentifier());
      for (const auto& source : run.getSources())
      {
        for (const auto& query : source.identifications)
        {
          if (!query.metaValueExists("old_map_index") && !query.metaValueExists("map_index")) continue;
          IdentificationData::Observation observation = query;
          if (observation.metaValueExists("old_map_index"))
          {
            Size old_map_index = (Size)observation.getMetaValue("old_map_index");
            Size file_index = (Size)observation.getMetaValue("map_index");
            observation.setMetaValue("map_index", mapid_table[make_pair(file_index, old_map_index)]);
            observation.removeMetaValue("old_map_index");
          }
          else
          {
            observation.removeMetaValue("map_index");
          }
          run.replaceObservation(query.getId(), observation);
        }
      }
    }
  }

  namespace
  {
    template <class MapType>
    FeatureGroupingAlgorithm::MapIdentifications mapIdentifications(const MapType& map)
    {
      FeatureGroupingAlgorithm::MapIdentifications result {map.getIdentificationData(), {}, {}};
      const auto collect = [&](const auto& self, const auto& feature) -> void {
        result.queries.insert(feature.getIDQueries().begin(), feature.getIDQueries().end());
        result.matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
        {
          for (const auto& subordinate : feature.getSubordinates())
          {
            self(self, subordinate);
          }
        }
      };
      for (const auto& feature : map)
      {
        collect(collect, feature);
      }
      return result;
    }

    template <class MapType>
    void groupMaps(const std::vector<MapType>& maps, ConsensusMap& grouped)
    {
      std::vector<FeatureGroupingAlgorithm::MapIdentifications> identifications;
      identifications.reserve(maps.size());
      for (const auto& map : maps)
      {
        identifications.push_back(FeatureGroupingAlgorithm::getMapIdentifications(map));
      }
      FeatureGroupingAlgorithm::groupIdentifications(std::move(identifications), grouped);
    }

    template <class MapType>
    void postprocessMaps(const std::vector<MapType>& maps, ConsensusMap& out)
    {
      if (std::any_of(maps.begin(), maps.end(), [](const MapType& map) { return IdentificationDataConverter::hasPeptideIdentifications(map); }))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Feature grouping needs the identifications of the maps as identification data; "
                                          "convert peptide identifications first (see groupWithIdentificationData_())");
      }
      FeatureGroupingAlgorithm::groupIdentifications(maps, out);

      // canonical ordering for checking the results:
      out.sortByQuality();
      out.sortByMaps();
      out.sortBySize();
    }

    template <class MapType>
    void groupWithData(const std::vector<MapType>& maps, ConsensusMap& out, const std::function<void(const std::vector<MapType>&)>& group)
    {
      const bool legacy = std::any_of(maps.begin(), maps.end(), [](const MapType& map) { return IdentificationDataConverter::hasPeptideIdentifications(map); });
      std::vector<MapType> converted;
      group(IdentificationDataConverter::withIdentificationData(maps, converted));
      if (legacy)
      {
        IdentificationDataConverter::exportConsensusIDs(out);
      }
    }
  } // namespace

  FeatureGroupingAlgorithm::MapIdentifications FeatureGroupingAlgorithm::getMapIdentifications(const FeatureMap& map)
  {
    return mapIdentifications(map);
  }

  FeatureGroupingAlgorithm::MapIdentifications FeatureGroupingAlgorithm::getMapIdentifications(const ConsensusMap& map)
  {
    return mapIdentifications(map);
  }

  void FeatureGroupingAlgorithm::groupIdentifications(std::vector<MapIdentifications> maps, ConsensusMap& grouped)
  {
    using ID = IdentificationData;
    std::set<ID::QueryReference> grouped_queries;
    std::set<ID::MatchReference> grouped_matches;
    for (const ConsensusFeature& feature : grouped)
    {
      grouped_queries.insert(feature.getIDQueries().begin(), feature.getIDQueries().end());
      grouped_matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
    }
    ID result;
    for (Size map_index = 0; map_index < maps.size(); ++map_index)
    {
      ID& data = maps[map_index].data;
      const auto& queries = maps[map_index].queries;
      const auto& matches = maps[map_index].matches;
      // Like the peptide identifications of features that grouping leaves out: what only such features link goes.
      const auto dropped = [&](const ID::MatchReference& match) {
        return matches.contains(match) && !grouped_matches.contains(match);
      };
      std::set<ID::QueryReference> drop;
      for (const auto& run : data.getRuns())
      {
        for (const auto& source : run.getSources())
        {
          for (const auto& query : source.identifications)
          {
            const ID::QueryReference reference {run.getUuid(), query.getId()};
            bool linked = queries.contains(reference), kept = grouped_queries.contains(reference), left = false;
            for (const auto& match : query.getMatches())
            {
              const ID::MatchReference match_reference {run.getUuid(), match.getId()};
              linked = linked || matches.contains(match_reference);
              kept = kept || grouped_matches.contains(match_reference);
              left = left || !dropped(match_reference);
            }
            if (linked && !kept && !left) drop.insert(reference);
          }
        }
      }
      data.eraseMatches([&](const ID::Run& run, const ID::Identification&, const ID::Match& match) {
        return dropped({run.getUuid(), match.getId()});
      });
      data.eraseIdentifications([&](const ID::Run& run, const ID::Identification& query) {
        return drop.contains({run.getUuid(), query.getId()});
      });
      for (const auto& current : data.getRuns())
      {
        auto& run = data.getRun(current.getIdentifier());
        for (const auto& source : run.getSources())
        {
          for (const auto& query : source.identifications)
          {
            ID::Observation observation = query;
            observation.setMetaValue("map_index", map_index);
            run.replaceObservation(query.getId(), observation);
          }
        }
      }
      result.merge(data);
    }
    grouped.getIdentificationData() = std::move(result);
  }

  void FeatureGroupingAlgorithm::groupIdentifications(const std::vector<FeatureMap>& maps, ConsensusMap& grouped)
  {
    groupMaps(maps, grouped);
  }

  void FeatureGroupingAlgorithm::groupIdentifications(const std::vector<ConsensusMap>& maps, ConsensusMap& grouped)
  {
    groupMaps(maps, grouped);
  }

  void FeatureGroupingAlgorithm::postprocess_(const std::vector<FeatureMap>& maps, ConsensusMap& out) const
  {
    postprocessMaps(maps, out);
  }

  void FeatureGroupingAlgorithm::postprocess_(const std::vector<ConsensusMap>& maps, ConsensusMap& out) const
  {
    postprocessMaps(maps, out);
  }

  void FeatureGroupingAlgorithm::groupWithIdentificationData_(const std::vector<FeatureMap>& maps, ConsensusMap& out,
                                                              const std::function<void(const std::vector<FeatureMap>&)>& group)
  {
    groupWithData(maps, out, group);
  }

  void FeatureGroupingAlgorithm::groupWithIdentificationData_(const std::vector<ConsensusMap>& maps, ConsensusMap& out,
                                                              const std::function<void(const std::vector<ConsensusMap>&)>& group)
  {
    groupWithData(maps, out, group);
  }

  FeatureGroupingAlgorithm::~FeatureGroupingAlgorithm() = default;

} //namespace OpenMS
