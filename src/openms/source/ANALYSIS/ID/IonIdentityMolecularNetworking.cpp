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
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <iostream>
#include <set>
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

namespace OpenMS
{
  void IonIdentityMolecularNetworking::annotateConsensusMap(ConsensusMap& out) const
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
    for (size_t i = 0; i < out.size(); i++)
    {
      out[i].setMetaValue(Constants::UserParam::IIMN_ROW_ID, i+1);
      if (!out[i].metaValueExists(Constants::UserParam::IIMN_LINKED_GROUPS)) continue;
      auto feature_vertex = add_vertex(VertexLabel(String(i), true), g);
      for (const auto& group: out[i].getMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS).toStringList())
      {
        if (!already_in_graph[group])
        {
          add_edge(feature_vertex, add_vertex(VertexLabel(group, false), g), g);
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
} // closing namespace OpenMS