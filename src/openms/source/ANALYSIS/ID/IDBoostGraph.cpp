// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/variant.hpp>
#include <boost/variant/detail/hash_variant.hpp>

using namespace OpenMS;
using namespace std;


namespace OpenMS
{
  void IDBoostGraph::buildGraph(const ProteinIdentification& proteins, const std::vector<PeptideIdentification>& psms)
  {
    //typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, IDBoostGraphNode> Graph;
    typedef boost::variant<const PeptideHit*, const ProteinHit*> IDPointer;
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, IDPointer> Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

    unordered_map<IDPointer, Graph::vertex_descriptor> vertex_map;

    unordered_set<string, ProteinHit*, ProteinHit::ProteinHitPtrAccessionHash> accession_set;

    vector<ProteinHit> protein_hits = proteins.getHits();
    std::transform (protein_hits.begin(), protein_hits.end(), accession_map.begin(), accession_map.begin(), std::plus<int>());

    Graph g;
    for (auto psm : psms)
    {
      //TODO add psm nodes here
      for (auto const& peptide : psm.getHits())
      {
        auto vertex_iter = vertex_map.find(&peptide);
        if (vertex_iter != vertex_map.end() )
        {
          g[pepnode] = iter->second;
        }
        else
        {
          vertex_t pepnode = boost::add_vertex(g);
          g[pepnode] = &peptide;
        }

        for (auto const& proteinAcc : peptide.extractProteinAccessionsSet())
        {

          for (auto const& protein : proteins.getHits())
          {
            vertex_t protnode = boost::add_vertex(g);
            g[protnode] = &protein;
          }
        }
      }
    }

  }
}

