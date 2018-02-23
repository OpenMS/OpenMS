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
#define INFERENCE_DEBUG

#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>

#include <ostream>


using namespace OpenMS;
using namespace std;


namespace OpenMS
{
  void IDBoostGraph::buildGraph_(const ProteinIdentification& proteins, const std::vector<PeptideIdentification>& psms)
  {
    //TODO add support for (consensus) feature information

    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};

    unordered_map<string, const ProteinHit*> accession_map{};

    std::transform (proteins.getHits().cbegin(), proteins.getHits().cend(), inserter(accession_map, accession_map.begin()),
                    [](const ProteinHit& p) {return make_pair<string, const ProteinHit*>(string(p.getAccession()), &p);});

    for (auto const& psm : psms)
    {
      //TODO add psm nodes here or take only the best hit!!
      for (auto const& peptide : psm.getHits())
      {
        vertex_t pepV = addVertexWithLookup_(&peptide, vertex_map);
        for (auto const& proteinAcc : peptide.extractProteinAccessionsSet())
        {
          auto const& accPHPair = accession_map.find(string(proteinAcc));
          vertex_t protV = addVertexWithLookup_(accPHPair->second, vertex_map);
          boost::add_edge(protV, pepV, g);
        }
      }
    }
  }

  /// Do sth on ccs
  void IDBoostGraph::doSomethingOnCC(const ProteinIdentification& protein, const std::vector<PeptideIdentification>& peptides){
    buildGraph_(protein, peptides);
    computeConnectedComponents_();
    for (unsigned int i = 0; i < numCCs; ++i)
    {
      boost::filtered_graph<Graph, boost::function<bool(edge_t)>, boost::function<bool(vertex_t)> > fg (g,
                                           [& i, this](edge_t e){return componentProperty[e.m_source] == i;},
                                           [& i, this](vertex_t v){return componentProperty[v] == i;});
      // TODO make function for writing graph?
      // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
      LabelVisitor lv;
      auto labels = boost::make_transform_value_property_map([lv](const IDPointer& p) {return boost::apply_visitor(lv, p);}, boost::get(boost::vertex_bundle, fg));

      std::cout << "Printing cc " << i << std::endl;
      //boost::print_graph(fg);
      boost::write_graphviz(std::cout, fg, boost::make_label_writer(labels));

      //------------------ Now actual inference ------------------- //
      // Skip cc without peptide or protein
      //TODO this currently does not work because we do not filter edges I think
      //TODO introduce edge types!
      //TODO or skip single nodes instead
      if (boost::num_edges(fg) >= 1)
      {
        //TODO make them parameters and/or estimate with gold/grid search
        MessagePasserFactory<unsigned long> mpf (0.9,0.01,0.5,1.0);
        BetheInferenceGraphBuilder<unsigned long> bigb;

        // Cluster peptides with same parents
        //TODO this could be sped up by a good hashing function for sets of uints and using unordered_map
        map< set<vertex_t>, set<vertex_t> > pepClusters; //maps the parent (protein) set to peptides that have the same
        boost::filtered_graph<Graph, boost::function<bool(edge_t)>, boost::function<bool(vertex_t)> >::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        for (; ui != ui_end; ++ui)
        {
          IDPointer curr_idObj = fg[*ui];
          if (curr_idObj.which() == 0) //it's a peptide
          {
            //TODO assert that there is at least one protein mapping to this peptide! Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            set<vertex_t> parents;
            boost::filtered_graph<Graph, boost::function<bool(edge_t)>, boost::function<bool(vertex_t)> >::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, fg);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (fg[*adjIt].which() == 1) //if there are only two types (pep,prot) this check for prot is unnecessary
              {
                parents.insert(*adjIt);
              }
            }

            auto clusterIt = pepClusters.find(parents);
            if (clusterIt != pepClusters.end())
            {
              clusterIt->second.insert(*ui);
            }
            else
            {
              pepClusters[parents] = set<vertex_t>({*ui});
            }
          }
          else if (curr_idObj.which() == 1)
          {
            //TODO allow prior probability here
            bigb.insert_dependency(mpf.createProteinFactor(*ui));
          }
        }

        int count = 1;
        for (auto const& setpair : pepClusters)
        {
          #ifdef INFERENCE_DEBUG
          for (auto const& j : setpair.first)
            std::cout << j << ",";

          std::cout << ": ";
          for (auto const& j : setpair.second)
            std::cout << j << ",";

          std::cout << std::endl;
          #endif

          unsigned long label = boost::num_vertices(fg) + (count++);
          bigb.insert_dependency(mpf.createPeptideProbabilisticAdderFactor(setpair.first, label));

          for (auto const& j : setpair.second) // foreach peptide
          {
            //TODO assert that the peptide score is of type PEP!
            bigb.insert_dependency(mpf.createSumEvidenceFactor(setpair.first.size(), label, j));
            IDPointer p = fg[j];
            bigb.insert_dependency(mpf.createPeptideEvidenceFactor(j, boost::get<const PeptideHit*>(p)->getScore()));
          }

        }
        InferenceGraph<unsigned long> ig = bigb.to_graph();
        //TODO setup scheduler, expose its params and start inference
        //TODO write results to the IDPointers (maybe I have to change them to non-const pointers?)
      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }

  }


  void IDBoostGraph::computeConnectedComponents_()
  {
    componentProperty.resize(num_vertices(g));
    numCCs = boost::connected_components(g, &componentProperty[0]);
  }

  IDBoostGraph::vertex_t IDBoostGraph::addVertexWithLookup_(const IDPointer& ptr, unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map)
  {
    vertex_t v;
    auto const& vertex_iter = vertex_map.find(ptr);
    if (vertex_iter != vertex_map.end() )
    {
      v = boost::vertex(vertex_iter->second, g);
    }
    else
    {
      v = boost::add_vertex(g);
      vertex_map[ptr] = v;
      g[v] = ptr;
    }
    return v;
  }
}

