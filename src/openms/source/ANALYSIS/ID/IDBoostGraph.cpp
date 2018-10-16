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
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>

#include <ostream>
#ifdef _OPENMP
#include <omp.h>
#endif

#define INFERENCE_DEBUG

using namespace OpenMS;
using namespace std;

//TODO go through the vectors and see if we can preallocate some.
namespace OpenMS
{
  struct MyUIntSetHasher
  {
  public:
    size_t operator()(const set<unsigned long>& s) const
        {
          return boost::hash_range(s.begin(), s.end());
        }
  };

  IDBoostGraph::PeptideCluster IDBoostGraph::staticPC{};
  IDBoostGraph::ProteinGroup IDBoostGraph::staticPG{};

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

    for (auto &prot : proteins_.getHits())
    {
      accession_map[prot.getAccession()] = &prot;
    }

    ProgressLogger pl;
    pl.setLogType(ProgressLogger::CMD);
    pl.startProgress(0,idedSpectra_.size(), "Building graph...");
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
          auto accToPHit = accession_map.find(std::string(proteinAcc));
          if (accToPHit == accession_map.end())
          {
            std::cout << "Warning: Building graph: skipping pep that maps to a non existent protein accession.";
            continue;
          }
          //TODO consider/calculate missing digests.
          //int missingTheorDigests = accToPHit->second->getMetaValue("missingTheorDigests");
          //accToPHit->second->setMetaValue("missingTheorDigests", missingTheorDigests);
          IDPointer prot(accToPHit->second);
          vertex_t protV = addVertexWithLookup_(prot, vertex_map);
          boost::add_edge(protV, pepV, g);
        }
      }
      pl.nextProgress();
    }
    pl.endProgress();
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
  void IDBoostGraph::applyFunctorOnCCs(std::function<void(Graph&)> functor)
  {
    if (ccs_.empty()) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }

    int i = 0;
    #pragma omp parallel for
    for (auto g_it = ccs_.begin() ; g_it < ccs_.end(); ++g_it)
    {
      #ifdef INFERENCE_DEBUG
      std::cout << omp_get_thread_num() << std::endl;
      std::cout << "Printing cc " << ++i << std::endl;
      printGraph(std::cout, *g_it);
      #endif

      functor(*g_it);
    }
  }

  void IDBoostGraph::annotateIndistProteins(bool addSingletons) const
  {
    if (ccs_.empty() && boost::num_vertices(g) == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Graph empty. Build it first.");
    }

    ProgressLogger pl;
    pl.setLogType(ProgressLogger::CMD);

    if (ccs_.empty())
    {
      pl.startProgress(0, 1, "Annotating indistinguishable proteins...");
      annotateIndistProteins_(g, addSingletons);
      pl.nextProgress();
      pl.endProgress();
    }
    else
    {
      int i = 0;

      #pragma omp parallel for
      for (auto g_it = ccs_.begin() ; g_it < ccs_.end(); ++g_it)
      {

        #ifdef INFERENCE_DEBUG
        std::cout << omp_get_thread_num() << std::endl;
        std::cout << "Printing cc " << i << std::endl;
        printGraph(std::cout, *g_it);
        #endif

        annotateIndistProteins_(*g_it, addSingletons);
        pl.setProgress(++i);
      }
      pl.endProgress();
    }
  }

  void IDBoostGraph::annotateIndistProteins_(const Graph& fg, bool addSingletons) const
  {
    //TODO evaluate hashing performance on sets
    unordered_map<PeptideNodeSet, ProteinNodeSet, MyUIntSetHasher > indistProteins; //find indist proteins

    Graph::vertex_iterator ui, ui_end;
    boost::tie(ui, ui_end) = boost::vertices(fg);

    // Cluster proteins
    for (; ui != ui_end; ++ui)
    {
      IDBoostGraph::IDPointer curr_idObj = fg[*ui];
      //TODO introduce an enum for the types to make it more clear.
      //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
      if (curr_idObj.which() == 0) //protein: find indist. ones
      {
        //TODO assert that there is at least one peptide mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
        //Or just check rigorously here.
        PeptideNodeSet childPeps;
        GraphConst::adjacency_iterator adjIt, adjIt_end;
        boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, fg);
        for (; adjIt != adjIt_end; ++adjIt)
        {
          if (fg[*adjIt].which() == 3) //if there are only two types (pep,prot) this check for pep is actually unnecessary
          {
            childPeps.insert(*adjIt);
          }
        }

        auto clusterIt = indistProteins.find(childPeps);
        if (clusterIt != indistProteins.end())
        {
          clusterIt->second.insert(*ui);
        }
        else
        {
          indistProteins[childPeps] = ProteinNodeSet({*ui});
        }
      }
    }

    // add the protein groups to the underlying ProteinGroup data structure only
    for (auto const &pepsToGrps : indistProteins)
    {
      if (pepsToGrps.second.size() <= 1 && !addSingletons)
      {
        continue;
      }

      ProteinIdentification::ProteinGroup pg{};

      pg.probability = -1.0;
      for (auto const &proteinVID : pepsToGrps.second)
      {
        ProteinHit *proteinPtr = boost::get<ProteinHit*>(fg[proteinVID]);
        pg.accessions.push_back(proteinPtr->getAccession());

        // the following sets the score of the group to the max
        // this might make not much sense if there was no inference yet -> score = 0
        // And one might also want to use other scoring systems
        // Anyway, without prior or add. info, all indist. proteins should have the same
        // score
        double oldscore = proteinPtr->getScore();
        if (oldscore > pg.probability)
        {
          pg.probability = oldscore;
        }
      }

      proteins_.getIndistinguishableProteins().push_back(pg);
    }
  }


  void IDBoostGraph::clusterIndistProteinsAndPeptides()
  {}
 /* void IDBoostGraph::clusterIndistProteinsAndPeptides()
  {
    if (numCCs_ == 0) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }

    #pragma omp parallel for
    for (unsigned int i = 0; i < numCCs_; ++i)
    {
      FilteredGraph fg(g,
                       [& i, this](edge_t e) { return componentProperty_[e.m_source] == i; },
                       [& i, this](vertex_t v) { return componentProperty_[v] == i; });

      #ifdef INFERENCE_DEBUG
      std::cout << "Printing cc " << i << std::endl;
      printFilteredGraph(std::cout, fg);
      #endif

      // Skip cc without peptide or protein
      //TODO this currently does not work because we do not filter edges I think
      //TODO introduce edge types or skip nodes without neighbors inside the if instead
      //TODO do quick bruteforce calculation if the cc is really small
      if (boost::num_vertices(fg) >= 3)
      {
        // Cluster peptides with same parents
        //TODO this could be sped up by a good hashing function for sets of uints and using unordered_map
        //TODO actually this version does not need the second set.
        unordered_map< ProteinNodeSet, PeptideNodeSet, MyUIntSetHasher > pepClusters; //maps the parent (protein) set to peptides that have the same
        unordered_map< PeptideNodeSet, ProteinNodeSet, MyUIntSetHasher > indistProteins; //find indist proteins

        FilteredGraph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        // Cluster proteins
        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::IDPointer curr_idObj = fg[*ui];
          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
          if (curr_idObj.which() == 0) //protein: find indist. ones
          {
            //TODO assert that there is at least one peptide mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            PeptideNodeSet childPeps;
            FilteredGraph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, fg);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (fg[*adjIt].which() == 3) //if there are only two types (pep,prot) this check for pep is actually unnecessary
              {
                childPeps.insert(*adjIt);
              }
            }

            auto clusterIt = indistProteins.find(childPeps);
            if (clusterIt != indistProteins.end())
            {
              clusterIt->second.insert(*ui);
            }
            else
            {
              indistProteins[childPeps] = ProteinNodeSet({*ui});
            }
          }
        }

        // add the protein groups to the graph
        // and edges from the groups to the proteins for quick access
        for (auto const& pepsToGrps : indistProteins)
        {
          if (pepsToGrps.second.size() <= 1)
            continue;
          //We can't point to protein groups while we fill them. Pointers invalidate in growing vectors.
          //proteins_.getIndistinguishableProteins().push_back(ProteinGroup{});
          //ProteinGroup& pg = proteins_.getIndistinguishableProteins().back();
          auto grpVID = boost::add_vertex(&staticPG, g);
          componentProperty_.push_back(i);
          for (auto const &proteinVID : pepsToGrps.second)
          {
            ProteinHit *proteinPtr = boost::get<ProteinHit*>(fg[proteinVID]);
            //pg.accessions.push_back(proteinPtr->getAccession());
            boost::add_edge(proteinVID, grpVID, g);
            for (auto const &pepVID : pepsToGrps.first)
            {
              boost::remove_edge(proteinVID, pepVID, g);
            }
          }
          for (auto const &pepVID : pepsToGrps.first)
          {
            boost::add_edge(grpVID, pepVID, g);
          }
          //pg.probability = -1.0;
        }


        // Create new filtering view and cluster peptides (TODO check if refiltering is actually necessary)
        // We won't come around resetting the iterators though, so not much saved.
        FilteredGraph fgNew(g,
                            [& i, this](edge_t e) { return componentProperty_[e.m_source] == i; },
                            [& i, this](vertex_t v) { return componentProperty_[v] == i; });
        FilteredGraph::vertex_iterator uiNew, uiNew_end;
        boost::tie(uiNew,uiNew_end) = boost::vertices(fgNew);

        for (; uiNew != uiNew_end; ++uiNew)
        {
          IDBoostGraph::IDPointer curr_idObj = fgNew[*uiNew];
          //TODO introduce an enum for the types to make it more clear.
          if (curr_idObj.which() == 3) //peptide: find peptide clusters
          {
            //TODO assert that there is at least one protein mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            ProteinNodeSet parents;
            IDBoostGraph::FilteredGraph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*uiNew, fgNew);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (fgNew[*adjIt].which() <= 1) // Either protein or protein group
              {
                parents.insert(*adjIt);
              }
            }

            auto clusterIt = pepClusters.find(parents);
            if (clusterIt != pepClusters.end())
            {
              clusterIt->second.insert(*uiNew);
            }
            else
            {
              pepClusters[parents] = PeptideNodeSet({*uiNew});
            }
          }
        }

        // we add an edge from protein to pepCluster and from pepCluster to peptides
        // peptides can use the same info from there.
        for (auto const& protsToPepClusters : pepClusters)
        {
          if (protsToPepClusters.first.size() <= 1)
            continue;
          auto pcVID = boost::add_vertex(&staticPC, g);
          componentProperty_.push_back(i);
          for (auto const& pgVID : protsToPepClusters.first)
          {
            boost::add_edge(pgVID, pcVID, g);
            for (auto const& peptideVID : protsToPepClusters.second)
            {
              boost::remove_edge(pgVID, peptideVID, g);
            }
          }
          for (auto const& peptideVID : protsToPepClusters.second)
          {
            boost::add_edge(pcVID, peptideVID, g);
          }
        }


        #ifdef INFERENCE_DEBUG
        std::cout << "Printing cc " << i << "with intermediate nodes." << std::endl;
        printFilteredGraph(std::cout, fgNew);
        #endif

      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  }*/


  //TODO we should probably rename it to splitCC now. Add logging and timing?
  void IDBoostGraph::computeConnectedComponents()
  {
    auto vis = dfs_ccsplit_visitor(ccs_);
    boost::depth_first_search(g, visitor(vis));
    LOG_INFO << "Found " << ccs_.size() << "CCs" << std::endl;
    g.clear();
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

  void IDBoostGraph::printFilteredGraph(std::ostream& out, const FilteredGraph& fg) const
  {
    // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
    //TODO Is passing "this" to lambda bad? How can I pass private members then?
    auto labels = boost::make_transform_value_property_map([this](const IDPointer &p) { return boost::apply_visitor(lv_, p); },
                                                           boost::get(boost::vertex_bundle, fg));
    //boost::print_graph(fg);
    boost::write_graphviz(out, fg, boost::make_label_writer(labels));
  }

  void IDBoostGraph::printGraph(std::ostream& out, const Graph& fg) const
  {
    // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
    //TODO Is passing "this" to lambda bad? How can I pass private members then?
    auto labels = boost::make_transform_value_property_map([this](const IDPointer &p) { return boost::apply_visitor(lv_, p); },
                                                           boost::get(boost::vertex_bundle, fg));
    //boost::print_graph(fg);
    boost::write_graphviz(out, fg, boost::make_label_writer(labels));
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

