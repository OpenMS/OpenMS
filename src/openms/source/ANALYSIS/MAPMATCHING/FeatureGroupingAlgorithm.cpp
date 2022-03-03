// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include <OpenMS/CONCEPT/Factory.h>

#include <unordered_set>
#include <regex>
#include <unordered_map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

// VertexLabel is used for IINM to build a bipartite graph using UndirectedOSMIdGraph
struct VertexLabel
{
    VertexLabel() = default;
    VertexLabel(OpenMS::String u, bool is_feat):uid(u), is_feature(is_feat) {}
    // if feature vertex: index (0 to n) for a feature in ConsensusMap; 
    // if group vertex: Constants::UserParam::ADDUCT_GROUP from ConsensusFeature Constants::UserParam::IIMN_LINKED_GROUPS
    OpenMS::String uid = "0"; 
    // false = vertex represents a group, true = vertex represents a feature   
    bool is_feature = false;
};

using UndirectedOSMIdGraph = boost::adjacency_list<
    boost::vecS,
    boost::vecS, 
    boost::undirectedS,
    VertexLabel>;

using namespace std;

namespace OpenMS
{
  //register products here
  void FeatureGroupingAlgorithm::registerChildren()
  {
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmLabeled::getProductName(), &FeatureGroupingAlgorithmLabeled::create);
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmUnlabeled::getProductName(), &FeatureGroupingAlgorithmUnlabeled::create);
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmQT::getProductName(), &FeatureGroupingAlgorithmQT::create);
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmKD::getProductName(), &FeatureGroupingAlgorithmKD::create);

  }

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
          Size old_map_index = id.getMetaValue("old_map_index");
          Size file_index = id.getMetaValue("map_index");
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
        Size old_map_index = id.getMetaValue("old_map_index");
        Size file_index = id.getMetaValue("map_index");
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

  void FeatureGroupingAlgorithm::annotateIonIdentityNetworks(ConsensusMap& out) const
  {
    // bipartite graph with ConsensusFeature indexes and Groups from Features
    // Vertexes contain uid (index/Group) and is_feature (bool)
    // e.g. 
    // ConsensusFeature0 contains Feature (Group = 1) and Feature (Group = 2)
    // ConsensusFeature1 contains Feature (Group = 2) and Feature (Group = 3)
    // ConsensusFeature2 contains Feature (Group = 4) and Feature (Group = 5)
    // graph looks like (Vertex of ConsensusFeature <--> Vertex of Group; component_number)
    // 0, true <--> 1, false; 0
    // 0, true <--> 2, false; 0
    // 1, true <--> 2, false; 0
    // 1, true <--> 3, false; 0
    // 2, true <--> 4, false; 1
    // 2, true <--> 5, false; 1
    // Total number of components: 2
    UndirectedOSMIdGraph g;

    // For each ConsensusFeature: add Constants::UserParam::IIMN_ROW_ID and if Constants::UserParam::IIMN_LINKED_GROUPS is present
    // add a Vertex for the ConsensusFeature and all corresponding Groups to bipartite graph.
    // Check if a Group Vertex has been added already, since the same Group can occur multiple times.
    for (size_t i = 0; i < out.size(); i++)
    {
      out[i].setMetaValue(Constants::UserParam::IIMN_ROW_ID, i+1);
      if (!out[i].metaValueExists(Constants::UserParam::IIMN_LINKED_GROUPS)) continue;
      auto feature_vertex = add_vertex(VertexLabel(String(i), true), g);
      for (const auto& group: out[i].getMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS).toStringList())
      {
        bool already_set = false;
        for (auto i : boost::make_iterator_range(vertices(g)))
        {
            if ((group == g[i].uid) && (!g[i].is_feature)) 
            {
              boost::add_edge(feature_vertex, i, g);
              already_set = true;
              break;
            }
        }
        if (!already_set)
        {
          boost::add_edge(feature_vertex, add_vertex(VertexLabel(group, false), g), g);
        }
      }
    }

    // represents the component number for each vertex
    std::vector<int> components (boost::num_vertices (g));
    boost::connected_components (g, &components[0]);

    // annotate network number and create a map with feature ID and partner IDs
    // partner feature vertexes are connected via a group vertex
    unordered_map<size_t, set<size_t>> partner_map;
    for (auto i : boost::make_iterator_range(vertices(g)))
    {
      if (!g[i].is_feature) continue;
      out[stoi(g[i].uid)].setMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER, components[i]+1);
      auto group_neighbours = boost::adjacent_vertices(i, g);
      for (auto gn : make_iterator_range(group_neighbours))
      {
        auto feature_partners = boost::adjacent_vertices(gn, g);
        for (auto partner : make_iterator_range(feature_partners))
        {
          if (i == partner) continue;
          partner_map[stoi(g[i].uid)].insert(stoi(g[partner].uid));
        }
      }
    }

    // annotate partners
    for (const auto& i : partner_map)
    {
      String partners;
      for (const auto& j : i.second)
      {
        if (partners.size() > 0) partners += ";";
        partners += out[j].getMetaValue(Constants::UserParam::IIMN_ROW_ID).toString();
      }
      out[i.first].setMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS, partners);
    }

    // remove Constants::UserParam::IIMN_LINKED_GROUPS meta values
    for (size_t i = 0; i < out.size(); i++)
    {
      out[i].removeMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS);
    }
 }

  FeatureGroupingAlgorithm::~FeatureGroupingAlgorithm()
  {
  }

} //namespace OpenMS
