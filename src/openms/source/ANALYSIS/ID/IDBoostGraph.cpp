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
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>

#include <ostream>
#ifdef _OPENMP
#include <omp.h>
#endif

//#define INFERENCE_DEBUG
#define INFERENCE_MT_DEBUG

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

  //TODO move to cpp file
  struct IDBoostGraph::SequenceToReplicateChargeVariantHierarchy
  {
    //TODO only add the intermediate nodes if there are more than one "splits"
    SequenceToReplicateChargeVariantHierarchy(Size nrReplicates, int minCharge, int maxCharge):
        seq_to_vecs_{},
        minCharge_(minCharge),
        nrCharges_(Size(maxCharge - minCharge) + 1u),
        nrReplicates_(nrReplicates)
    {}

    void insert(String& seq, Size replicate, int charge, vertex_t pepVtx)
    {
      auto seq_it = seq_to_vecs_.emplace(std::move(seq), std::vector<std::vector<std::set<vertex_t>>>{nrReplicates_, {nrCharges_, std::set<vertex_t>()}});
      seq_it.first->second[replicate][charge - minCharge_].insert(pepVtx);
    }

    void insertToGraph(vertex_t rootProteinVtx, Graph& graph)
    {

      for (const auto& seqContainer : seq_to_vecs_)
      {
        vertex_t pep = boost::add_vertex(Peptide{seqContainer.first}, graph);

        vector<vertex_t> prots_for_pepseq;
        GraphConst::adjacency_iterator adjIt, adjIt_end;
        boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(rootProteinVtx, graph);
        for (; adjIt != adjIt_end; adjIt++)
        {
          IDBoostGraph::IDPointer curr_idObj = graph[*adjIt];
          if (curr_idObj.which() == 0) //protein
          {
            boost::add_edge(*adjIt, pep, graph);
            prots_for_pepseq.push_back(*adjIt);
          }
        }

        for (Size s = 0; s < seqContainer.second.size(); ++s)
        {
          vertex_t ri = boost::add_vertex(RunIndex{s},graph);
          boost::add_edge(pep, ri, graph);
          for (Size t = 0; t < seqContainer.second[s].size(); ++t)
          {
            vertex_t cs = boost::add_vertex(Charge{minCharge_ + int(t)}, graph);
            boost::add_edge(ri, cs, graph);
            for (const auto& pepVtx : seqContainer.second[s][t])
            {
              // TODO: Instead of collecting and here in the innermost loop removing
              // the edges from proteins to PSMs, we could (if we assume that nothing else
              // was added yet) just clear all (in-)edges for the pepVtx here
              for (const auto& parent : prots_for_pepseq)
              {
                boost::remove_edge(parent, pepVtx, graph);
              }

              boost::add_edge(cs, pepVtx, graph);
            }
          }
        }
      }
    }

    std::unordered_map<std::string, std::vector<std::vector<std::set<vertex_t>>>> seq_to_vecs_;

    int minCharge_;
    Size nrCharges_;
    Size nrReplicates_;
  };

  IDBoostGraph::IDBoostGraph(ProteinIdentification& proteins,
                             std::vector<PeptideIdentification>& idedSpectra):
      proteins_(proteins),
      idedSpectra_(idedSpectra),
      exp_design_(ExperimentalDesign::fromIdentifications({proteins}))
  {}

  IDBoostGraph::IDBoostGraph(ProteinIdentification& proteins,
      std::vector<PeptideIdentification>& idedSpectra,
      const ExperimentalDesign& ed):
    proteins_(proteins),
    idedSpectra_(idedSpectra),
    exp_design_(ed)
  {}

  //TODO actually to build the graph, the inputs could be passed const. But if you want to do sth
  //on the graph later it needs to be non-const. Overload this function or somehow make sure it can be used const.
  void IDBoostGraph::buildGraph(Size use_top_psms, bool readstore_run_info)
  {
    StringList runs;
    proteins_.getPrimaryMSRunPath(runs);
    Size nrRuns = runs.size();
    //TODO add support for (consensus) feature information

    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};

    unordered_map<string, ProteinHit*> accession_map{};

    for (auto& prot : proteins_.getHits())
    {
      accession_map[prot.getAccession()] = &prot;
    }

    ProgressLogger pl;
    pl.setLogType(ProgressLogger::CMD);
    pl.startProgress(0, idedSpectra_.size(), "Building graph...");
    for (auto& spectrum : idedSpectra_)
    {
      Size run(0);
      //TODO check if the spectrum is from one of the runs in the ProtID? or assert that?
      if (readstore_run_info)
      {
        if (spectrum.metaValueExists("map_index"))
        {
          run = spectrum.getMetaValue("map_index");
          if (run >= nrRuns)
          {
            std::cout << "Warning: Reference to non-existing run found at peptide ID. Skipping." << std::endl;
            continue;
          }
        }
        else
        {
          std::cout << "Warning: Trying to read run information but none present at peptide ID. Skipping." << std::endl;
          continue;
        }
      }
      //TODO add psm regularizer nodes here optionally if using multiple psms
      auto pepIt = spectrum.getHits().begin();
      //TODO sort or assume sorted
      auto pepItEnd = use_top_psms == 0 || spectrum.getHits().empty() ? spectrum.getHits().end() : spectrum.getHits().begin() + use_top_psms;
      for (; pepIt != pepItEnd; ++pepIt)
      {
        IDPointer pepPtr(&(*pepIt));
        vertex_t pepV = addVertexWithLookup_(pepPtr, vertex_map);
        pepHitVtx_to_run_[pepV] = run;

        for (auto const & proteinAcc : pepIt->extractProteinAccessionsSet())
        {
          // assumes protein is present
          auto accToPHit = accession_map.find(std::string(proteinAcc));
          if (accToPHit == accession_map.end())
          {
            std::cout << "Warning: Building graph: skipping pep that maps to a non existent protein accession.";
            continue;
          }
          //TODO consider/calculate missing digests. Probably not here though!
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

  //TODO actually to build the graph, the inputs could be passed const. But if you want to do sth
  //on the graph later it needs to be non-const. Overload this function or somehow make sure it can be used const.
  /*
  template <class T>
  void IDBoostGraph::buildTheoreticalGraph(pair<int,int> chargeRange, unsigned int nrReplicates, FASTAContainer<T>& proteins)
  {
    ///ProteaseDigestion enzyme; //TODO get from protein ID run
    bool IL_equivalent = true; //TODO how to incorporate that?
    // cache the first proteins
    const size_t PROTEIN_CACHE_SIZE = 4e5; // 400k should be enough for most DB's and is not too hard on memory either (~200 MB FASTA)

    proteins.cacheChunk(PROTEIN_CACHE_SIZE);

    if (proteins.empty()) // we do not allow an empty database
    {
      LOG_ERROR << "Error: An empty database was provided. Mapping makes no sense. Aborting..." << std::endl;
      //TODO throw Exception
    }

    bool has_active_data = true;
    bool invalid_protein_sequence = false;
    Size count_j_proteins(0);
    Size prot_count(0);
    String prot = "";
    while (true)
    {
      has_active_data = proteins.activateCache(); // swap in last cache

      if (!has_active_data)
        break; // leave while-loop

      proteins.cacheChunk(PROTEIN_CACHE_SIZE);


      for (Size i = 0; i < prot_count; ++i)
      {
        prot = proteins.chunkAt(i).sequence;
        prot.remove('*');

        // check for invalid sequences with modifications
        if (prot.has('[') || prot.has('('))
        {
          invalid_protein_sequence = true; // not omp-critical because its write-only
          // we cannot throw an exception here, since we'd need to catch it within the parallel region
        }

        // convert  L/J to I; also replace 'J' in proteins
        if (IL_equivalent)
        {
          prot.substitute('L', 'I');
          prot.substitute('J', 'I');
        }
        else
        { // warn if 'J' is found (it eats into aaa_max)
          if (prot.has('J'))
          {
            ++count_j_proteins;
          }
        }

        Size prot_idx = i + proteins.getChunkOffset();
      }

    }
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
  }*/

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

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
    {
      Graph& curr_cc = ccs_.at(i);
      LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;

      #ifdef INFERENCE_MT_DEBUG
      LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
      #endif

      #ifdef INFERENCE_DEBUG
      LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, *g_it);
      LOG_INFO << "Printed cc " << i << std::endl;
      #endif

      functor(curr_cc);
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
      pl.startProgress(0, ccs_.size(), "Annotating indistinguishable proteins...");
      #pragma omp parallel for
      for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
      {
        const Graph& curr_cc = ccs_.at(i);
        LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;

        #ifdef INFERENCE_MT_DEBUG
        LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
        #endif

        #ifdef INFERENCE_DEBUG
        LOG_INFO << "Printing cc " << i << std::endl;
        printGraph(LOG_INFO, *g_it);
        LOG_INFO << "Printed cc " << i << std::endl;
        #endif

        annotateIndistProteins_(curr_cc, addSingletons);
        pl.setProgress(i);
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
          if (fg[*adjIt].which() >= 3) //if there are only two types (pep,prot) this check for pep is actually unnecessary
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

      // TODO you could allocate as many groups as proteins in the beginning
      // then you do not need a critical section. Resize afterwards.
      #pragma omp critical (ProteinGroups)
      {proteins_.getIndistinguishableProteins().push_back(pg);};
    }
  }


  /*void IDBoostGraph::clusterIndistProteinsAndPeptidesOld()
  {
    if (ccs_.empty()) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
    {
      Graph& curr_cc = ccs_.at(i);

      LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;

      #ifdef INFERENCE_MT_DEBUG
      LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
      #endif

      #ifdef INFERENCE_DEBUG
      LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, curr_cc);
      LOG_INFO << "Printed cc " << i << std::endl;
      #endif

      // Skip cc without peptide or protein
      //TODO better to do quick bruteforce calculation if the cc is really small
      if (boost::num_edges(curr_cc) >= 1)
      {
        // Cluster peptides with same parents
        unordered_map< ProteinNodeSet, PeptideNodeSet, MyUIntSetHasher > pepClusters; //maps the parent (protein) set to peptides that have the same
        unordered_map< PeptideNodeSet, ProteinNodeSet, MyUIntSetHasher > indistProteins; //find indist proteins

        Graph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(curr_cc);

        // Cluster proteins
        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::IDPointer curr_idObj = curr_cc[*ui];

          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
          if (curr_idObj.which() == 0) //protein: find indist. ones
          {
            //TODO assert that there is at least one peptide mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            PeptideNodeSet childPeps;
            Graph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, curr_cc);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (curr_cc[*adjIt].which() == 3) //if there are only two types (pep,prot) this check for pep is actually unnecessary
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
          auto grpVID = boost::add_vertex(&staticPG, curr_cc);

          for (auto const &proteinVID : pepsToGrps.second)
          {
            //ProteinHit *proteinPtr = boost::get<ProteinHit*>(curr_cc[proteinVID]);
            //pg.accessions.push_back(proteinPtr->getAccession());
            boost::add_edge(proteinVID, grpVID, curr_cc);
            for (auto const &pepVID : pepsToGrps.first)
            {
              boost::remove_edge(proteinVID, pepVID, curr_cc);
            }
          }
          for (auto const &pepVID : pepsToGrps.first)
          {
            boost::add_edge(grpVID, pepVID, curr_cc);
          }
          //pg.probability = -1.0;
        }



        // reset iterator to loop through vertices again for peptide clusters
        boost::tie(ui,ui_end) = boost::vertices(curr_cc);

        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::IDPointer curr_idObj = curr_cc[*ui];
          //TODO introduce an enum for the types to make it more clear.
          if (curr_idObj.which() == 3) //peptide: find peptide clusters
          {
            //TODO assert that there is at least one protein mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            ProteinNodeSet parents;
            Graph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, curr_cc);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (curr_cc[*adjIt].which() <= 1) // Either protein or protein group
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
              pepClusters[parents] = PeptideNodeSet({*ui});
            }
          }
        }

        // we add an edge from protein to pepCluster and from pepCluster to peptides
        // peptides can use the same info from there.
        for (auto const& protsToPepClusters : pepClusters)
        {
          if (protsToPepClusters.first.size() <= 1)
            continue;
          auto pcVID = boost::add_vertex(&staticPC, curr_cc);
          for (auto const& pgVID : protsToPepClusters.first)
          {
            boost::add_edge(pgVID, pcVID, curr_cc);
            for (auto const& peptideVID : protsToPepClusters.second)
            {
              boost::remove_edge(pgVID, peptideVID, curr_cc);
            }
          }
          for (auto const& peptideVID : protsToPepClusters.second)
          {
            boost::add_edge(pcVID, peptideVID, curr_cc);
          }
        }

        #ifdef INFERENCE_DEBUG
        LOG_INFO << "Printing cc " << i << "with intermediate nodes." << std::endl;
        printGraph(LOG_INFO, curr_cc);
        LOG_INFO << "Printed cc " << i << "with intermediate nodes." << std::endl;
        #endif

      }
      else
      {
        LOG_INFO << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  }*/

  void IDBoostGraph::clusterIndistProteinsAndPeptides()
  {
    Size nrReplicates = 1;
    if (!pepHitVtx_to_run_.empty()) //graph built with run info
    {
      StringList runs;
      proteins_.getPrimaryMSRunPath(runs);
      nrReplicates = runs.size();
    }

    pair<int,int> chargeRange = proteins_.getSearchParameters().getChargeRange();

    if (ccs_.empty()) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
    {
      Graph& curr_cc = ccs_[i];

      LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;

      #ifdef INFERENCE_MT_DEBUG
      LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
      #endif

      #ifdef INFERENCE_DEBUG
      LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, curr_cc);
      LOG_INFO << "Printed cc " << i << std::endl;
      #endif

      // Skip cc without peptide or protein
      //TODO better to do quick bruteforce calculation if the cc is really small
      if (boost::num_edges(curr_cc) >= 1)
      {

        Graph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(curr_cc);

        // Cluster peptides with same sequence and create a replicate and charge hierarchy underneath
        for (; ui != ui_end; ++ui)
        {
          SequenceToReplicateChargeVariantHierarchy hierarchy{nrReplicates, chargeRange.first, chargeRange.second};
          if (curr_cc[*ui].which() == 0) //protein: same seq peptideHits have to be at a single protein
          {
            Graph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, curr_cc);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              //pepHit, this also makes sure that pepHits already in hierarchy are masked
              if (curr_cc[*adjIt].which() == 3)
              {
                PeptideHit *phitp = boost::get<PeptideHit *>(curr_cc[*adjIt]);
                String seq = phitp->getSequence().toUnmodifiedString();

                Size rep = 0; //In case no replicate info was read.
                if (!pepHitVtx_to_run_.empty()) rep = pepHitVtx_to_run_[*adjIt];
                int chg = phitp->getCharge();

                hierarchy.insert(seq, rep, chg, *adjIt);
              }
            }
            hierarchy.insertToGraph(*ui, g);
          }
        }

        /*boost::tie(ui,ui_end) = boost::vertices(curr_cc);

         //TODO FINISH OR MOVE TO BUILD_THEORETICAL_GRAPH
        // Cluster peptides with same sequence
        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::IDPointer curr_idObj = curr_cc[*ui];
          unordered_map<string, map<pair<unsigned int, int>, set<vertex_t>>> seq_to_rep_chg_to_pepVariantNodes{};
          if (curr_idObj.which() == 0) //protein: same seq peptideHits have to be at a single protein
          {
            Graph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, curr_cc);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (curr_cc[*adjIt].which() == 3) //pepHit
              {
                PeptideHit* phitp = boost::get<PeptideHit*>(curr_cc[*adjIt]);
                String seq = phitp->getSequence().toUnmodifiedString();
                //TODO we could specify/check in the beginning if it is a multirun experiment and how many runs it has
                //TODO make sure (prob. best during merging) that different map_indices indeed represent replicates
                //or pass a lookup table here to get from map_index to replicate.
                int rep = phitp->metaValueExists("map_index") ? int(phitp->getMetaValue("map_index")) : 0;
                int chg = phitp->getCharge();
                auto seq_it = seq_to_rep_chg_to_pepVariantNodes.find(seq);
                if (seq_it == seq_to_rep_chg_to_pepVariantNodes.end())
                {
                  map<pair<unsigned int, int>, set<vertex_t>> initMap{};
                  initMap.insert({{rep,chg}, set<vertex_t>{*adjIt}});
                  seq_to_rep_chg_to_pepVariantNodes.insert({seq, initMap});
                }
                else
                {
                  auto rep_chg_it = seq_it->second.find({rep,chg});
                  if (rep_chg_it == seq_it->second.end())
                  {
                    seq_it->second.insert({{rep,chg}, set<vertex_t>{*adjIt}});
                  }
                  else
                  {
                    rep_chg_it->second.insert(*adjIt);
                  }
                }
              }
            }
            //TODO finish converting the hierarchy into graph nodes.
            //TODO only add the intermediate nodes if there are more than one "splits"
            for (auto seq_to_map : seq_to_rep_chg_to_pepVariantNodes)
            {
              auto seqVID = boost::add_vertex(seq_to_map.first, curr_cc);
              for (int poss_rep : poss_reps)
              {
                auto begin = seq_to_map.second.lower_bound({poss_rep, INT_MIN});
                auto end = seq_to_map.second.upper_bound({poss_rep, INT_MAX});
                if (begin != end) // we have more than one replicate
                {

                }
                auto repVID = boost::add_vertex(seq_to_map.first, curr_cc);
              }
              for (auto rep_chg_to_vtcs : seq_to_map.second)
              {
                //TODO remove the edge from ui (prot) to the peptides in the set
              }
            }

          }
        }*/

        // Cluster peptides with same parents
        unordered_map< ProteinNodeSet, PeptideNodeSet, MyUIntSetHasher > pepClusters; //maps the parent (protein) set to peptides that have the same
        unordered_map< PeptideNodeSet, ProteinNodeSet, MyUIntSetHasher > indistProteins; //find indist proteins

        boost::tie(ui,ui_end) = boost::vertices(curr_cc);

        // Cluster proteins
        for (; ui != ui_end; ++ui)
        {
          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
          if (curr_cc[*ui].which() == 0) //protein: find indist. ones
          {
            //TODO assert that there is at least one peptide mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            PeptideNodeSet childPeps;
            Graph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, curr_cc);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (curr_cc[*adjIt].which() == 3) //if there are only two types (pep,prot) this check for pep is actually unnecessary
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
          auto grpVID = boost::add_vertex(ProteinGroup{}, curr_cc);

          for (auto const &proteinVID : pepsToGrps.second)
          {
            //ProteinHit *proteinPtr = boost::get<ProteinHit*>(curr_cc[proteinVID]);
            //pg.accessions.push_back(proteinPtr->getAccession());
            boost::add_edge(proteinVID, grpVID, curr_cc);
            for (auto const &pepVID : pepsToGrps.first)
            {
              boost::remove_edge(proteinVID, pepVID, curr_cc);
            }
          }
          for (auto const &pepVID : pepsToGrps.first)
          {
            boost::add_edge(grpVID, pepVID, curr_cc);
          }
          //pg.probability = -1.0;
        }



        // reset iterator to loop through vertices again for peptide clusters
        boost::tie(ui,ui_end) = boost::vertices(curr_cc);

        for (; ui != ui_end; ++ui)
        {
          //TODO introduce an enum for the types to make it more clear.
          if (curr_cc[*ui].which() == 3) //peptide: find peptide clusters
          {
            //TODO assert that there is at least one protein mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            ProteinNodeSet parents;
            Graph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, curr_cc);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (curr_cc[*adjIt].which() <= 1) // Either protein or protein group
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
              pepClusters[parents] = PeptideNodeSet({*ui});
            }
          }
        }

        // we add an edge from protein to pepCluster and from pepCluster to peptides
        // peptides can use the same info from there.
        for (auto const& protsToPepClusters : pepClusters)
        {
          if (protsToPepClusters.first.size() <= 1)
            continue;
          auto pcVID = boost::add_vertex(PeptideCluster{}, curr_cc);
          for (auto const& pgVID : protsToPepClusters.first)
          {
            boost::add_edge(pgVID, pcVID, curr_cc);
            for (auto const& peptideVID : protsToPepClusters.second)
            {
              boost::remove_edge(pgVID, peptideVID, curr_cc);
            }
          }
          for (auto const& peptideVID : protsToPepClusters.second)
          {
            boost::add_edge(pcVID, peptideVID, curr_cc);
          }
        }

        #ifdef INFERENCE_DEBUG
        LOG_INFO << "Printing cc " << i << "with intermediate nodes." << std::endl;
        printGraph(LOG_INFO, curr_cc);
        LOG_INFO << "Printed cc " << i << "with intermediate nodes." << std::endl;
        #endif

      }
      else
      {
        LOG_INFO << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  }


  //TODO we should probably rename it to splitCC now. Add logging and timing?
  void IDBoostGraph::computeConnectedComponents()
  {
    auto vis = dfs_ccsplit_visitor(ccs_);
    boost::depth_first_search(g, visitor(vis));
    LOG_INFO << "Found " << ccs_.size() << " connected components." << std::endl;
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

  ////// Hashers for the strong typedefs
  std::size_t hash_value(const IDBoostGraph::Peptide& x)
  {
    boost::hash<std::string> hasher;
    return hasher(static_cast<std::string>(x));
  }
  std::size_t hash_value(const IDBoostGraph::RunIndex& x)
  {
    boost::hash<Size> hasher;
    return hasher(static_cast<Size>(x));
  }
  std::size_t hash_value(const IDBoostGraph::Charge& x)
  {
    boost::hash<int> hasher;
    return hasher(static_cast<int>(x));
  }
  std::size_t hash_value(const IDBoostGraph::ProteinGroup&)
  {
    return 0;
  }
  std::size_t hash_value(const IDBoostGraph::PeptideCluster&)
  {
    return 1;
  }

}

