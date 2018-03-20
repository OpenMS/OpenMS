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
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>

#include <ostream>


using namespace OpenMS;
using namespace std;

//TODO go through the vectors and see if we can preallocate some.
namespace OpenMS
{

  IDBoostGraph::IDBoostGraph(ProteinIdentification& proteins, std::vector<PeptideIdentification>& idedSpectra):
    proteins_(proteins),
    idedSpectra_(idedSpectra)
  {}


  //TODO actually to build the graph, the inputs could be passed const. But if you want to do sth
  //on the graph later it needs to be non-const. Overload this function or somehow make sure it can be used const.
  void IDBoostGraph::buildGraph(bool use_all_psms)
  {
    //TODO add support for (consensus) feature information

    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};

    unordered_map<string, ProteinHit*> accession_map{};

    auto protIt = proteins_.getHits().begin();
    auto protItEnd = proteins_.getHits().end();
    for (; protIt != protItEnd; ++protIt)
    {
      accession_map[string(protIt->getAccession())] = &(*protIt);
    }

    for (auto & spectrum : idedSpectra_)
    {
      //TODO add psm nodes here if using all psms
      auto pepIt = spectrum.getHits().begin();
      auto pepItEnd = use_all_psms || spectrum.getHits().empty() ? spectrum.getHits().end() : spectrum.getHits().begin() + 1;
      for (; pepIt != pepItEnd; ++pepIt)
      {
        IDPointer pepPtr(&(*pepIt));
        vertex_t pepV = addVertexWithLookup_(pepPtr, vertex_map);
        for (auto const & proteinAcc : pepIt->extractProteinAccessionsSet())
        {
          // assumes protein is present
          IDPointer prot(accession_map.find(std::string(proteinAcc))->second);
          vertex_t protV = addVertexWithLookup_(prot, vertex_map);
          boost::add_edge(protV, pepV, g);
        }
      }
    }
  }

/* Const version
 * void IDBoostGraph::buildGraph(const ProteinIdentification& proteins, const std::vector<PeptideIdentification>& psms)
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
  void IDBoostGraph::applyFunctorOnCCs(std::function<void(FilteredGraph&)> functor)
  {
    if (numCCs_ == 0) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }
    for (unsigned int i = 0; i < numCCs_; ++i)
    {
      FilteredGraph fg(g,
                       [& i, this](edge_t e) { return componentProperty_[e.m_source] == i; },
                       [& i, this](vertex_t v) { return componentProperty_[v] == i; });

      #ifdef INFERENCE_DEBUG
      // TODO make function for writing graph?
      // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
      LabelVisitor lv;
      auto labels = boost::make_transform_value_property_map([lv](IDPointer &p) { return boost::apply_visitor(lv, p); },
                                                             boost::get(boost::vertex_bundle, fg));
      std::cout << "Printing cc " << i << std::endl;
      //boost::print_graph(fg);
      boost::write_graphviz(std::cout, fg, boost::make_label_writer(labels));
      #endif

      functor(fg);
    }
  }

  /// Do sth on ccs
  void IDBoostGraph::annotateIndistinguishableGroups()
  {
    if (numCCs_ == 0) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }
    for (unsigned int i = 0; i < numCCs_; ++i)
    {
      FilteredGraph fg(g,
                       [& i, this](edge_t e) { return componentProperty_[e.m_source] == i; },
                       [& i, this](vertex_t v) { return componentProperty_[v] == i; });

      #ifdef INFERENCE_DEBUG
      // TODO make function for writing graph?
      // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
      LabelVisitor lv;
      auto labels = boost::make_transform_value_property_map([lv](IDPointer &p) { return boost::apply_visitor(lv, p); },
                                                             boost::get(boost::vertex_bundle, fg));
      std::cout << "Printing cc " << i << std::endl;
      //boost::print_graph(fg);
      boost::write_graphviz(std::cout, fg, boost::make_label_writer(labels));
      #endif

      // Skip cc without peptide or protein
      //TODO this currently does not work because we do not filter edges I think
      //TODO introduce edge types or skip nodes without neighbors inside the if instead
      //TODO do quick bruteforce calculation if the cc is really small
      if (boost::num_vertices(fg) >= 3)
      {
        // Cluster peptides with same parents
        //TODO this could be sped up by a good hashing function for sets of uints and using unordered_map
        map< set<IDBoostGraph::vertex_t>, set<IDBoostGraph::vertex_t> > pepClusters; //maps the parent (protein) set to peptides that have the same
        boost::filtered_graph<IDBoostGraph::Graph, boost::function<bool(IDBoostGraph::edge_t)>, boost::function<bool(IDBoostGraph::vertex_t)> >::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);


        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::IDPointer curr_idObj = fg[*ui];
          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
          if (curr_idObj.which() == 0) //it's a peptide
          {
            //TODO assert that there is at least one protein mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            set<IDBoostGraph::vertex_t> parents;
            IDBoostGraph::FilteredGraph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, fg);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (fg[*adjIt].which() == 1) //if there are only two types (pep,prot) this check for prot is actually unnecessary
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
              pepClusters[parents] = std::set<IDBoostGraph::vertex_t>({*ui});
            }
          }
        }

        proteins_.getIndistinguishableProteins().clear();
        // Afterwards go through all pepClusters, i.e. indist. groups and add them
        for (auto const& setpair : pepClusters)
        {
          ProteinIdentification::ProteinGroup pg;
          auto grp = boost::add_vertex(&pg, g);
          for (auto const& proteinVID : setpair.first)
          {
            ProteinHit* proteinPtr = boost::get<ProteinHit*>(fg[proteinVID]);
            pg.accessions.push_back(proteinPtr->getAccession());
            boost::add_edge(grp, proteinVID, g);
          }
          pg.probability = -1.0;
          proteins_.getIndistinguishableProteins().push_back(pg);
        }

      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  }


  void IDBoostGraph::computeConnectedComponents()
  {
    componentProperty_.resize(num_vertices(g));
    numCCs_ = boost::connected_components(g, &componentProperty_[0]);
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

/* Const version
 * IDBoostGraph::vertex_t IDBoostGraph::addVertexWithLookup_(IDPointerConst& ptr, unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>>& vertex_map)
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

