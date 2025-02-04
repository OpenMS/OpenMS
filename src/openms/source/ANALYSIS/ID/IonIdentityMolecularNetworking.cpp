// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

// VertexLabel is used for IINM to build a bipartite graph using UndirectedOSMIdGraph
struct VertexLabel
{
  VertexLabel() = default;
  VertexLabel(size_t u, bool is_feat):uid(u), is_feature(is_feat) {}
  // if feature vertex: index (0 to n) for a feature in ConsensusMap;
  // if group vertex: Constants::UserParam::ADDUCT_GROUP from ConsensusFeature Constants::UserParam::IIMN_LINKED_GROUPS
  size_t uid = 0; 
  // false = vertex represents a group, true = vertex represents a feature   
  bool is_feature = false;
};

using UndirectedOSMIdGraph = boost::adjacency_list<
      boost::vecS,
      boost::vecS, 
      boost::undirectedS,
      VertexLabel>;

namespace OpenMS
{
  void IonIdentityMolecularNetworking::annotateConsensusMap(ConsensusMap& consensus_map)
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
    // Check if a Group Vertex has been added already to the graph, since the same Group can occur multiple times.
    std::unordered_map<String, size_t> already_in_graph; // <group_uid, vertex_index>
    for (size_t i = 0; i < consensus_map.size(); i++)
    {
      consensus_map[i].setMetaValue(Constants::UserParam::IIMN_ROW_ID, i+1);
      if (!consensus_map[i].metaValueExists(Constants::UserParam::IIMN_LINKED_GROUPS)) continue;
      auto feature_vertex = add_vertex(VertexLabel(i, true), g);
      for (const auto& group: consensus_map[i].getMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS).toStringList())
      {
        if (!already_in_graph[group])
        {
          add_edge(feature_vertex, add_vertex(VertexLabel(std::stoull(group), false), g), g);
          already_in_graph[group] = boost::num_vertices(g);
          continue;
        }
        else add_edge(feature_vertex, already_in_graph[group]-1, g);
      }
    }

    // represents the component number for each vertex
    std::vector<int> components (boost::num_vertices(g));
    boost::connected_components (g, &components[0]);

    // annotate network number and create a map with feature ID and partner IDs
    // partner feature vertexes are connected via a group vertex
    std::unordered_map<size_t, std::set<size_t>> partner_map;
    for (const auto& i : boost::make_iterator_range(vertices(g)))
    {
      if (!g[i].is_feature) continue;
      consensus_map[g[i].uid].setMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER, components[i]+1);
      auto group_neighbours = boost::adjacent_vertices(i, g);
      for (const auto& gn : make_iterator_range(group_neighbours))
      {
        auto feature_partners = boost::adjacent_vertices(gn, g);
        
        for (const auto& partner : make_iterator_range(feature_partners))
        {
          if (i == partner) continue;
          partner_map[g[i].uid].insert(g[partner].uid);
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
        partners += consensus_map[j].getMetaValue(Constants::UserParam::IIMN_ROW_ID).toString();
      }
      consensus_map[i.first].setMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS, partners);
    }

    // remove Constants::UserParam::IIMN_LINKED_GROUPS meta values
    for (size_t i = 0; i < consensus_map.size(); i++)
    {
      consensus_map[i].removeMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS);
    }
  }
  /**
  @brief Generates a supplementary pairs table required for GNPS IIMN, as defined here: https://ccms-ucsd.github.io/GNPSDocumentation/fbmn-iin/#supplementary-pairs
  */
  void IonIdentityMolecularNetworking::writeSupplementaryPairTable(const ConsensusMap& consensus_map, const String& output_file)
  {
    // exit early if there is no IIMN annotations (first feature has no Constants::UserParam::IIMN_ROW_ID)
    if (!consensus_map[0].metaValueExists(Constants::UserParam::IIMN_ROW_ID)) return;

    // generate unordered map with feature id and partner feature ids
    std::unordered_map<size_t, std::set<size_t>> feature_partners; // map<feature_index, partner_feature_indices>
    for (size_t i = 0; i < consensus_map.size(); i++)
    {
      if (!consensus_map[i].metaValueExists(Constants::UserParam::IIMN_ADDUCT_PARTNERS)) continue;
      std::stringstream ss(consensus_map[i].getMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS).toChar());
      while(ss.good())
      {
        String substr;
        getline(ss, substr, ';');
        feature_partners[i].insert(std::stoi(substr)-1);
      }
    }

    // sort feature partners in new map and swap values (this will give reproducible results across different runs)
    std::map<size_t, std::set<size_t>> sorted;
    for (auto& [key, value] : feature_partners)
    {
    sorted[key].swap(value);
    }

    // get number of partners for each feature to later calculate score of annotation
    std::unordered_map<size_t, size_t> num_partners;
    for (const auto& entry: sorted)
    {
      num_partners[entry.first] = entry.second.size();
    }

    // initialize SVOutStream with tab separation
    std::ofstream outstr(output_file.c_str());
    SVOutStream out(outstr, ",", "_", String::NONE);
    
    // write table header
    out << "ID1" << "ID2" << "EdgeType" << "Score" << "Annotation" << std::endl;

    // write edge annotation for each feature / partner feature pair (sorted)
    for (const auto& entry: sorted)
    {
      for (const auto& partner_index: entry.second)
      {
        out << consensus_map[entry.first].getMetaValue(Constants::UserParam::IIMN_ROW_ID);
        out << consensus_map[partner_index].getMetaValue(Constants::UserParam::IIMN_ROW_ID);
        out << "MS1 annotation";
        out << num_partners[entry.first] + num_partners[partner_index] - 2; // total number of direct partners from both features minus themselves
        std::stringstream annotation;
        annotation << consensus_map[entry.first].getMetaValue(Constants::UserParam::IIMN_BEST_ION, String("default")) << " "
                   << consensus_map[partner_index].getMetaValue(Constants::UserParam::IIMN_BEST_ION, String("default")) << " "
                   << "dm/z=" << String(std::abs(consensus_map[entry.first].getMZ() - consensus_map[partner_index].getMZ()));
        out << annotation.str();
        out << std::endl;
        sorted[partner_index].erase(entry.first); // remove other direction to avoid duplicates
      }
    }
    outstr.close();
  }
} // closing namespace OpenMS