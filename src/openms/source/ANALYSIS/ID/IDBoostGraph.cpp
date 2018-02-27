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
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/bind.hpp>
#include <boost/graph/connected_components.hpp>

#include <ostream>


using namespace OpenMS;
using namespace std;

//TODO go through the vectors and see if we can preallocate some.
namespace OpenMS
{
  //TODO actually to build the graph, the inputs could be passed const. But if you want to do sth
  //on the graph later it needs to be non-const. Overload this function or somehow make sure it can be used const.
  void IDBoostGraph::buildGraph_(ProteinIdentification& proteins, std::vector<PeptideIdentification>& psms)
  {
    //TODO add support for (consensus) feature information

    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};

    unordered_map<string, ProteinHit*> accession_map{};

    auto protIt = proteins.getHits().begin();
    auto protItEnd = proteins.getHits().end();
    for (; protIt != protItEnd; ++protIt)
    {
      accession_map[string(protIt->getAccession())] = &(*protIt);
    }

    for (auto & psm : psms)
    {
      //TODO add psm nodes here or take only the best hit!!
      auto psmIt = psm.getHits().begin();
      auto psmItEnd = psm.getHits().end();
      for (; psmIt != psmItEnd; ++psmIt)
      {
        IDPointer pep(&(*psmIt));
        vertex_t pepV = addVertexWithLookup_(pep, vertex_map);
        for (auto const & proteinAcc : psmIt->extractProteinAccessionsSet())
        {
          // assumes protein is present
          IDPointer prot(accession_map.find(std::string(proteinAcc))->second);
          vertex_t protV = addVertexWithLookup_(prot, vertex_map);
          boost::add_edge(protV, pepV, g);
        }
      }
    }
  }

/*  void IDBoostGraph::buildGraph_(const ProteinIdentification& proteins, const std::vector<PeptideIdentification>& psms)
  {
    //TODO add support for (consensus) feature information

    unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>> vertex_map{};

    unordered_map<string, const ProteinHit*> accession_map{};

    // choose cbegin and cend if you do const
    std::transform (proteins.getHits().cbegin(), proteins.getHits().cend(), inserter(accession_map, accession_map.begin()),
                    [](const ProteinHit& p) {return make_pair<string, const ProteinHit*>(string(p.getAccession()), &p);});

    // add const after auto for const version
    for (auto const & psm : psms)
    {
      //TODO add psm nodes here or take only the best hit!!
      for (auto const & peptide : psm.getHits())
      {
        IDPointerConst pep = &peptide;
        vertex_t pepV = addVertexWithLookup_(pep, vertex_map);
        for (auto const & proteinAcc : peptide.extractProteinAccessionsSet())
        {
          // assumes protein is present
          IDPointerConst prot = accession_map.find(std::string(proteinAcc))->second;
          vertex_t protV = addVertexWithLookup_(prot, vertex_map);
          boost::add_edge(protV, pepV, gconst);
        }
      }
    }
  }*/

  /// Do sth on ccs
  void IDBoostGraph::applyFunctorOnCCs(ProteinIdentification &protein,
                                       std::vector<PeptideIdentification> &peptides,
                                       std::function<void(FilteredGraph &)> functor)
  {
    buildGraph_(protein, peptides);
    computeConnectedComponents_();
    for (unsigned int i = 0; i < numCCs; ++i)
    {
      FilteredGraph fg(g,
                       [& i, this](edge_t e) { return componentProperty[e.m_source] == i; },
                       [& i, this](vertex_t v) { return componentProperty[v] == i; });
      // TODO make function for writing graph?
      // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
      LabelVisitor lv;
      auto labels = boost::make_transform_value_property_map([lv](IDPointer &p) { return boost::apply_visitor(lv, p); },
                                                             boost::get(boost::vertex_bundle, fg));

      std::cout << "Printing cc " << i << std::endl;
      //boost::print_graph(fg);
      boost::write_graphviz(std::cout, fg, boost::make_label_writer(labels));

      functor(fg);
    }
  }

/*  void IDBoostGraph::applyFunctorOnCCs(ProteinIdentification& protein, std::vector<PeptideIdentification>& peptides){
    buildGraph_(protein, peptides);
    computeConnectedComponents_();
    for (unsigned int i = 0; i < numCCs; ++i)
    {
      FilteredGraph fg (g,
                        [& i, this](edge_t e){return componentProperty[e.m_source] == i;},
                        [& i, this](vertex_t v){return componentProperty[v] == i;});
      // TODO make function for writing graph?
      // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
      LabelVisitor lv;
      auto labels = boost::make_transform_value_property_map([lv](IDPointer& p) {return boost::apply_visitor(lv, p);}, boost::get(boost::vertex_bundle, fg));

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

        // Store the IDs of the nodes for which you want the posteriors in the end (usually at least proteins)
        // Maybe later peptides (e.g. for an iterative procedure)
        vector<vector<unsigned long>> posteriorVars;

        for (; ui != ui_end; ++ui)
        {
          IDPointer curr_idObj = fg[*ui];
          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
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
            //TODO allow an already present prior probability here
            bigb.insert_dependency(mpf.createProteinFactor(*ui));
            posteriorVars.push_back({*ui});
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
            bigb.insert_dependency(mpf.createPeptideEvidenceFactor(j, boost::get<PeptideHit*>(p)->getScore()));
          }

        }
        InferenceGraph<unsigned long> ig = bigb.to_graph();
        //TODO we should parametrize the type of scheduler and its params.
        PriorityScheduler<unsigned long> scheduler(0.001, 1e-8, 1ul<<32);
        scheduler.add_ab_initio_edges(ig);
        BeliefPropagationInferenceEngine<unsigned long> bpie(scheduler, ig);
        auto posteriorFactors = bpie.estimate_posteriors(posteriorVars);

        for (auto const& posteriorFactor : posteriorFactors)
        {
          double posterior = 0.0;
          SetPosteriorVisitor pv;
          unsigned long nodeId = posteriorFactor.ordered_variables()[0];
          const PMF& pmf = posteriorFactor.pmf();
          // If Index 1 is in the range of this result PMFFactor it is non-zero
          if (1 >= pmf.first_support()[0] && 1 <= pmf.last_support()[0]) {
            posterior = pmf.table()[1 - pmf.first_support()[0]];
          }
          auto bound_visitor = std::bind(pv, std::placeholders::_1, posterior);
          boost::apply_visitor(bound_visitor, fg[nodeId]);
        }
      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }

  }*/


  void IDBoostGraph::computeConnectedComponents_()
  {
    componentProperty.resize(num_vertices(g));
    numCCs = boost::connected_components(g, &componentProperty[0]);
  }

  IDBoostGraph::vertex_t IDBoostGraph::addVertexWithLookup_(IDPointer& ptr, unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map)
  {
    vertex_t v;
    auto vertex_iter = vertex_map.find(ptr);
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

/*  IDBoostGraph::vertex_t IDBoostGraph::addVertexWithLookup_(IDPointerConst& ptr, unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>>& vertex_map)
  {
    vertex_t v;
    auto const& vertex_iter = vertex_map.find(ptr);
    if (vertex_iter != vertex_map.end() )
    {
      v = boost::vertex(vertex_iter->second, gconst);
    }
    else
    {
      v = boost::add_vertex(gconst);
      vertex_map[ptr] = v;
      gconst[v] = ptr;
    }
    return v;
  }*/
}

