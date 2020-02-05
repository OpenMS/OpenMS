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
#include <OpenMS/SYSTEM/StopWatch.h>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>

#include <ostream>
#ifdef _OPENMP
#include <omp.h>
#endif

//#define INFERENCE_DEBUG

//#define INFERENCE_MT_DEBUG

using namespace OpenMS;
using namespace std;
using Internal::IDBoostGraph;

//TODO go through the vectors and see if we can preallocate some.
namespace OpenMS
{

  /// Hasher for sets of uints using boost::hash_range
  struct MyUIntSetHasher
  {
  public:
    size_t operator()(const set<IDBoostGraph::vertex_t>& s) const
        {
          return boost::hash_range(s.begin(), s.end());
        }
  };

  /// Helper struct to create Sequence->Replicate->Chargestate hierarchy for a set of PSMs from a protein
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
      int chargeToPut = charge - minCharge_;
      OPENMS_PRECONDITION(replicate < nrReplicates_, "Replicate OOR")
      OPENMS_PRECONDITION(static_cast<Size>(chargeToPut) < nrCharges_, "Charge OOR")

      auto seq_it = seq_to_vecs_.emplace(std::move(seq), std::vector<std::vector<std::set<vertex_t>>>{nrReplicates_, std::vector<std::set<vertex_t>>(nrCharges_, std::set<vertex_t>())});
      seq_it.first->second[replicate][chargeToPut].insert(pepVtx);
    }

    //TODO finish and rework (root not needed?)
    void insertToGraph(vertex_t /*rootProteinVtx*/, Graph& graph)
    {
      for (const auto& seqContainer : seq_to_vecs_)
      {
        vertex_t pep = boost::add_vertex(Peptide{seqContainer.first}, graph);

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
              GraphConst::adjacency_iterator adjIt, adjIt_end;
              // This assumes, that at this point, only proteins are connected to PSMS
              boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(pepVtx, graph);
              for (; adjIt != adjIt_end; adjIt++)
              {
                IDBoostGraph::IDPointer curr_idObj = graph[*adjIt];
                if (curr_idObj.which() == 0) //protein
                {
                  //below would invalidate iterator. We use clear vertex below
                  //boost::remove_edge(*adjIt, pepVtx, graph); //remove old one from protein
                  boost::add_edge(*adjIt, pep, graph); //instead add it to the sequence
                }
              }
              boost::clear_vertex(pepVtx, graph);
              boost::add_edge(cs, pepVtx, graph); //now connect the last level (charges) to this spectrum
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
                             std::vector<PeptideIdentification>& idedSpectra,
                             Size use_top_psms,
                             bool use_run_info,
                             const boost::optional<const ExperimentalDesign>& ed):
      protIDs_(proteins)
  {
    OPENMS_LOG_INFO << "Building graph on " << idedSpectra.size() << " spectra and " << proteins.getHits().size() << " proteins." << std::endl;
    if (use_run_info)
    {
      buildGraphWithRunInfo_(proteins, idedSpectra, use_top_psms, ed.get_value_or(ExperimentalDesign::fromIdentifications({proteins})));
    }
    else
    {
      buildGraph_(proteins, idedSpectra, use_top_psms);
    }
  }

  IDBoostGraph::IDBoostGraph(ProteinIdentification& proteins,
                             ConsensusMap& cmap,
                             Size use_top_psms,
                             bool use_run_info,
                             bool use_unassigned_ids,
                             const boost::optional<const ExperimentalDesign>& ed):
      protIDs_(proteins)
  {
    OPENMS_LOG_INFO << "Building graph on " << cmap.size() << " features, " << cmap.getUnassignedPeptideIdentifications().size() <<
    " unassigned spectra (if chosen) and " << proteins.getHits().size() << " proteins." << std::endl;
    if (use_run_info)
    {
      buildGraphWithRunInfo_(proteins, cmap, use_top_psms, use_unassigned_ids, ed.get_value_or(ExperimentalDesign::fromConsensusMap(cmap)));
    }
    else
    {
      buildGraph_(proteins, cmap, use_top_psms, use_unassigned_ids);
    }
  }

  unordered_map<unsigned, unsigned> convertMapLabelFree_(
      const map<pair<String, unsigned>, unsigned>& fileToRun,
      const StringList& files)
  {
    unordered_map<unsigned, unsigned> indexToRun;
    unsigned i = 0;
    for (const auto& file : files)
    {
      indexToRun[i] = fileToRun.at({file,1});
      ++i;
    } // TODO what if file is not in the experimental design? Check in the very beginning!?
    return indexToRun;
  }

  unordered_map<unsigned, unsigned> convertMap_(
      const map<pair<String, unsigned>, unsigned>& fileLabToPrefractionationGroup,
      const ConsensusMap::ColumnHeaders& idxToFileLabMappings,
      const String& experiment_type)
  {
    unordered_map<unsigned, unsigned> indexToRun;
    for (const auto& mapping : idxToFileLabMappings)
    {
      indexToRun[mapping.first] =
          fileLabToPrefractionationGroup.at(make_pair(mapping.second.filename, mapping.second.getLabelAsUInt(experiment_type)));
    } // TODO what if file is not in the experimental design? Check in the very beginning!?
    return indexToRun;
  }


  void IDBoostGraph::addPeptideIDWithAssociatedProteins_(PeptideIdentification& spectrum,
      unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map,
      const unordered_map<string, ProteinHit*>& accession_map,
      Size use_top_psms)
  {
    //TODO add psm regularizer nodes here optionally if using multiple psms (i.e. forcing them, so that only 1 or maybe 2 are present per spectrum)
    auto pepIt = spectrum.getHits().begin();
    //TODO sort or assume sorted
    auto pepItEnd = (use_top_psms == 0 || (spectrum.getHits().size() <= use_top_psms)) ? spectrum.getHits().end() : spectrum.getHits().begin() + use_top_psms;
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
         OPENMS_LOG_WARN << "Warning: Building graph: skipping pep that maps to a non existent protein accession." << std::endl;
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
  }

  void IDBoostGraph::addPeptideAndAssociatedProteinsWithRunInfo_(
      PeptideIdentification& spectrum,
      unordered_map<unsigned, unsigned>& indexToPrefractionationGroup,
      unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map,
      unordered_map<std::string, ProteinHit*>& accession_map,
      Size use_top_psms)
  {
    Size idx(0);
    Size pfg(0);

    if (spectrum.metaValueExists("map_index"))
    {
      idx = spectrum.getMetaValue("map_index");
      auto find_it = indexToPrefractionationGroup.find(idx);
      if (find_it == indexToPrefractionationGroup.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Reference (map_index) to non-existing run found at peptide ID."
            " Sth went wrong during merging. Aborting.");
      }
      pfg = find_it->second - 1; // Experimental design numbering starts at one
    }
    else
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Trying to read run information (map_index) but none present at peptide ID."
        " Did you annotate runs during merging? Aborting.");
    }

    //TODO add psm regularizer nodes here optionally if using multiple psms
    auto pepIt = spectrum.getHits().begin();
    //TODO sort or assume sorted
    auto pepItEnd = use_top_psms == 0 || spectrum.getHits().empty() ? spectrum.getHits().end() : spectrum.getHits().begin() + use_top_psms;
    for (; pepIt != pepItEnd; ++pepIt)
    {
      IDPointer pepPtr(&(*pepIt));
      vertex_t pepV = addVertexWithLookup_(pepPtr, vertex_map);

      //------- Only difference to the function without run info -----//
      pepHitVtx_to_run_[pepV] = pfg;
      //------- Only difference to the function without run info -----//

      for (auto const & proteinAcc : pepIt->extractProteinAccessionsSet())
      {
        // assumes protein is present
        auto accToPHit = accession_map.find(std::string(proteinAcc));
        if (accToPHit == accession_map.end())
        {
         OPENMS_LOG_WARN << "Warning: Building graph: skipping pep that maps to a non existent protein accession." << std::endl;
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
  }

  void IDBoostGraph::buildGraphWithRunInfo_(ProteinIdentification& proteins,
                                           ConsensusMap& cmap,
                                           Size use_top_psms,
                                           bool use_unassigned_ids,
                                           const ExperimentalDesign& ed)
  {
    unordered_map<unsigned, unsigned> indexToPrefractionationGroup;

    {
      // TODO check that the files in the ProteinID run are all in the Exp. Design
      //StringList files;
      //proteins.getPrimaryMSRunPath(files); // files merged in the protein identification run to be inferred
      const ConsensusMap::ColumnHeaders& colHeaders = cmap.getColumnHeaders(); // all possible files and labels in the experiment
      //TODO use exp. design to merge fractions
      map<pair<String, unsigned>, unsigned> fileLabelToPrefractionationGroup = ed.getPathLabelToPrefractionationMapping(false);
      nrPrefractionationGroups_ = fileLabelToPrefractionationGroup.size();
      indexToPrefractionationGroup = convertMap_(fileLabelToPrefractionationGroup, colHeaders, cmap.getExperimentType()); // convert to index in the peptide ids
    }

    //TODO is this vertex_map really necessary. I think PSMs are always unique in our datastructures and could be
    // added without lookup.
    // And for the proteins we could add the vertex ID to the accession_map here and use that for lookup
    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};
    unordered_map<std::string, ProteinHit*> accession_map{};

    for (auto& prot : proteins.getHits())
    {
      accession_map[prot.getAccession()] = &prot;
    }

    ProgressLogger pl;
    Size roughNrIds = cmap.size();
    if (use_unassigned_ids) roughNrIds += cmap.getUnassignedPeptideIdentifications().size();
    pl.setLogType(ProgressLogger::CMD);
    pl.startProgress(0, roughNrIds, "Building graph with run information...");
    const String& protRun = proteins.getIdentifier();
    for (auto& feat : cmap)
    {
      for (auto& spectrum : feat.getPeptideIdentifications())
      {
        if (spectrum.getIdentifier() == protRun)
        {
          addPeptideAndAssociatedProteinsWithRunInfo_(spectrum, indexToPrefractionationGroup,
                                                      vertex_map, accession_map, use_top_psms);
        }
      }
      pl.nextProgress();
    }

    if (use_unassigned_ids)
    {
      for (auto& id : cmap.getUnassignedPeptideIdentifications())
      {
        if (id.getIdentifier() == protRun)
        {
          addPeptideAndAssociatedProteinsWithRunInfo_(id, indexToPrefractionationGroup,
                                                      vertex_map, accession_map, use_top_psms);
        }
        pl.nextProgress();
      }
    }
    pl.endProgress();
  }

  void IDBoostGraph::buildGraphWithRunInfo_(ProteinIdentification& proteins,
                                           std::vector<PeptideIdentification>& idedSpectra,
                                           Size use_top_psms,
                                           const ExperimentalDesign& ed)
  {
    unordered_map<unsigned, unsigned> indexToPrefractionationGroup;

    {
      StringList files;
      proteins.getPrimaryMSRunPath(files);
      map<pair<String, unsigned>, unsigned> fileLabelToPrefractionationGroup = ed.getPathLabelToPrefractionationMapping(false);
      nrPrefractionationGroups_ = fileLabelToPrefractionationGroup.size();
      //TODO if only given proteins and peptide IDs we automatically assume label-free since I don't know
      // where the label would be stored.
      indexToPrefractionationGroup = convertMapLabelFree_(fileLabelToPrefractionationGroup, files); // convert to index in the peptide ids
    }

    //TODO is this vertex_map really necessary. I think PSMs are always unique in our datastructures and could be
    // added without lookup.
    // And for the proteins we could add the vertex ID to the accession_map here and use that for lookup
    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};
    unordered_map<std::string, ProteinHit*> accession_map{};

    for (auto& prot : proteins.getHits())
    {
      accession_map[prot.getAccession()] = &prot;
    }

    ProgressLogger pl;
    pl.setLogType(ProgressLogger::CMD);
    pl.startProgress(0, idedSpectra.size(), "Building graph with run info...");
    const String& protRun = proteins.getIdentifier();
    for (auto& spectrum : idedSpectra)
    {
      if (spectrum.getIdentifier() == protRun)
      {
        addPeptideAndAssociatedProteinsWithRunInfo_(spectrum, indexToPrefractionationGroup,
                                                      vertex_map, accession_map, use_top_psms);
      }
      pl.nextProgress();
    }
    pl.endProgress();
  }

  //TODO actually to build the graph, the inputs could be passed const. But if you want to do sth
  // on the graph later it needs to be non-const. Overload the next functions or somehow make sure it can be used const.
  void IDBoostGraph::buildGraph_(ProteinIdentification& proteins,
                                std::vector<PeptideIdentification>& idedSpectra,
                                Size use_top_psms)
  {
    StringList runs;
    proteins.getPrimaryMSRunPath(runs);

    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};

    unordered_map<string, ProteinHit*> accession_map{};

    for (auto& prot : proteins.getHits())
    {
      accession_map[prot.getAccession()] = &prot;
    }

    ProgressLogger pl;
    pl.setLogType(ProgressLogger::CMD);
    pl.startProgress(0, idedSpectra.size(), "Building graph...");
    const String& protRun = proteins.getIdentifier();
    for (auto& spectrum : idedSpectra)
    {
      if (spectrum.getIdentifier() == protRun)
      {
        addPeptideIDWithAssociatedProteins_(spectrum, vertex_map, accession_map, use_top_psms);
      }
      pl.nextProgress();
    }
    pl.endProgress();
  }


  void IDBoostGraph::buildGraph_(ProteinIdentification& proteins,
                                 ConsensusMap& cmap,
                                 Size use_top_psms,
                                 bool use_unassigned_ids)
  {
    StringList runs;
    proteins.getPrimaryMSRunPath(runs);

    unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>> vertex_map{};

    unordered_map<string, ProteinHit*> accession_map{};

    for (auto& prot : proteins.getHits())
    {
      accession_map[prot.getAccession()] = &prot;
    }

    ProgressLogger pl;
    Size roughNrIds = cmap.size();
    if (use_unassigned_ids) roughNrIds += cmap.getUnassignedPeptideIdentifications().size();
    pl.setLogType(ProgressLogger::CMD);
    pl.startProgress(0, roughNrIds, "Building graph...");
    const String& protRun = proteins.getIdentifier();
    for (auto& feature : cmap)
    {
      for (auto& id : feature.getPeptideIdentifications())
      {
        if (id.getIdentifier() == protRun)
        {
          addPeptideIDWithAssociatedProteins_(id, vertex_map, accession_map, use_top_psms);
        }
      }
      pl.nextProgress();
    }
    if (use_unassigned_ids)
    {
      for (auto& id : cmap.getUnassignedPeptideIdentifications())
      {
        if (id.getIdentifier() == protRun)
        {
          addPeptideIDWithAssociatedProteins_(id, vertex_map, accession_map, use_top_psms);
        }
        pl.nextProgress();
      }
    }
    pl.endProgress();
  }


  /* This would be a version where you try to build the graph based on the theoretical peptides
   * but this is quite some work and additional memory overhead
   * For now my plan is to create theoretically indistinguishable groups and (maybe nr of only?) missing peptides
   * in PeptideIndexer and use them here

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
     OPENMS_LOG_ERROR << "Error: An empty database was provided. Mapping makes no sense. Aborting..." << std::endl;
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


  /// Do sth on ccs
  void IDBoostGraph::applyFunctorOnCCs(const std::function<unsigned long(Graph&)>& functor)
  {
    if (ccs_.empty()) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }

    // Use dynamic schedule because big CCs take much longer!
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
    {
      #ifdef INFERENCE_BENCH
      StopWatch sw;
      sw.start();
      #endif

      Graph& curr_cc = ccs_.at(i);

      #ifdef INFERENCE_MT_DEBUG
     OPENMS_LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
      #endif

      #ifdef INFERENCE_DEBUG
     OPENMS_LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;
     OPENMS_LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, curr_cc);
     OPENMS_LOG_INFO << "Printed cc " << i << std::endl;
      #endif

      #ifdef INFERENCE_BENCH
      unsigned long result = functor(curr_cc);
      #else
      functor(curr_cc);
      #endif


      #ifdef INFERENCE_BENCH
      sw.stop();
      sizes_and_times_[i] = tuple<vertex_t, vertex_t, unsigned long, double>{boost::num_vertices(curr_cc), boost::num_edges(curr_cc), result, sw.getClockTime()};
      #endif
    }

    #ifdef INFERENCE_BENCH
    ofstream debugfile;
    debugfile.open("idgraph_functortimes_" + DateTime::now().getTime() + ".tsv");

    for (const auto& size_time : sizes_and_times_ )
    {
      debugfile << std::get<0>(size_time) << "\t" << std::get<1>(size_time) <<  "\t" << std::get<2>(size_time) << "\t" << std::get<3>(size_time) << "\n";
    }
    debugfile.close();
    #endif
  }

  /// Do sth on ccs single-threaded
  void IDBoostGraph::applyFunctorOnCCsST(const std::function<void(Graph&)>& functor)
  {
    if (ccs_.empty()) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }

    for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
    {
      #ifdef INFERENCE_BENCH
      StopWatch sw;
      sw.start();
      #endif

      Graph& curr_cc = ccs_.at(i);

      #ifdef INFERENCE_DEBUG
     OPENMS_LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;
     OPENMS_LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, curr_cc);
     OPENMS_LOG_INFO << "Printed cc " << i << std::endl;
      #endif

      functor(curr_cc);

      #ifdef INFERENCE_BENCH
      sw.stop();
      sizes_and_times_[i] = tuple<vertex_t, vertex_t, unsigned long, double>{boost::num_vertices(curr_cc), boost::num_edges(curr_cc), 0, sw.getClockTime()};
      #endif
    }

    #ifdef INFERENCE_BENCH
    ofstream debugfile;
    debugfile.open("idgraph_functortimes_" + DateTime::now().getTime() + ".tsv");

    for (const auto& size_time : sizes_and_times_ )
    {
      debugfile << std::get<0>(size_time) << "\t" << std::get<1>(size_time) <<  "\t" << std::get<2>(size_time) << "\t" << std::get<3>(size_time) << "\n";
    }
    debugfile.close();
    #endif
  }

  void IDBoostGraph::annotateIndistProteins(bool addSingletons)
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

        #ifdef INFERENCE_MT_DEBUG
       OPENMS_LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
        #endif

        #ifdef INFERENCE_DEBUG
       OPENMS_LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;
       OPENMS_LOG_INFO << "Printing cc " << i << std::endl;
        printGraph(LOG_INFO, curr_cc);
       OPENMS_LOG_INFO << "Printed cc " << i << std::endl;
        #endif

        annotateIndistProteins_(curr_cc, addSingletons);
        pl.setProgress(i);
      }
      pl.endProgress();
    }
  }

  void IDBoostGraph::calculateAndAnnotateIndistProteins(bool addSingletons)
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

        #ifdef INFERENCE_MT_DEBUG
       OPENMS_LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
        #endif

        #ifdef INFERENCE_DEBUG
       OPENMS_LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;
       OPENMS_LOG_INFO << "Printing cc " << i << std::endl;
        printGraph(LOG_INFO, curr_cc);
       OPENMS_LOG_INFO << "Printed cc " << i << std::endl;
        #endif

        calculateAndAnnotateIndistProteins_(curr_cc, addSingletons);
        pl.setProgress(i);
      }
      pl.endProgress();
    }
  }

  void IDBoostGraph::calculateAndAnnotateIndistProteins_(const Graph& fg, bool addSingletons)
  {
    //TODO evaluate hashing performance on sets
    unordered_map<PeptideNodeSet, ProteinNodeSet, MyUIntSetHasher> indistProteins; //find indist proteins

    Graph::vertex_iterator ui, ui_end;
    boost::tie(ui, ui_end) = boost::vertices(fg);

    //TODO refactor into function
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
      //  then you do not need a critical section. Resize afterwards.
      //  Or make a local vector of Groups and merge in a single threaded section
      #pragma omp critical (ProteinGroups)
      {protIDs_.getIndistinguishableProteins().push_back(pg);}
    }
  }

  void IDBoostGraph::annotateIndistProteins_(const Graph& fg, bool addSingletons)
  {
    Graph::vertex_iterator ui, ui_end;
    boost::tie(ui,ui_end) = boost::vertices(fg);

    for (; ui != ui_end; ++ui)
    {
      if (fg[*ui].which() == 1) //prot group
      {
        ProteinIdentification::ProteinGroup pg{};
        pg.probability = (double) boost::get<IDBoostGraph::ProteinGroup>(fg[*ui]); //init
        Graph::adjacency_iterator nbIt, nbIt_end;
        boost::tie(nbIt, nbIt_end) = boost::adjacent_vertices(*ui, fg);

        ProteinHit *proteinPtr = nullptr;
        for (; nbIt != nbIt_end; ++nbIt)
        {
          if (fg[*nbIt].which() == 0) //neighboring proteins
          {
            proteinPtr = boost::get<ProteinHit*>(fg[*nbIt]);
            pg.accessions.push_back(proteinPtr->getAccession());
          }
        }
        if (addSingletons || pg.accessions.size() > 1)
        {
          // TODO you could allocate as many groups as proteins in the beginning
          //  then you do not need a critical section. Resize afterwards.
          //  Or make a local vector of Groups and merge in a single threaded section
          #pragma omp critical (ProteinGroups)
          {protIDs_.getIndistinguishableProteins().push_back(pg);}
        }
      }
    }
  }

  void IDBoostGraph::getUpstreamNodesNonRecursive(std::queue<vertex_t>& q, Graph graph, int lvl, bool stop_at_first, std::vector<vertex_t>& result)
  {
    while (!q.empty())
    {
      vertex_t curr_node = q.front();
      q.pop();
      Graph::adjacency_iterator adjIt, adjIt_end;
      boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(curr_node, graph);
      for (;adjIt != adjIt_end; ++adjIt)
      {
        if (graph[*adjIt].which() <= lvl)
        {
          result.emplace_back(*adjIt);
          if (!stop_at_first && graph[*adjIt].which() > graph[curr_node].which())
          {
            q.emplace(*adjIt);
          }
        }
        else if (graph[*adjIt].which() > graph[curr_node].which())
        {
            q.emplace(*adjIt);
        }
      }
    }
  }

  /* Under development
  void IDBoostGraph::resolveGraphProteinCentric_(const Graph& fg)
  {
    Graph::vertex_iterator ui, ui_end;
    boost::tie(ui,ui_end) = boost::vertices(fg);

    unordered_set<vertex_t> prots;
    for (; ui != ui_end; ++ui)
    {
      if (fg[*ui].which() == 0) //prot
      {
        Graph::adjacency_iterator nbIt, nbIt_end;
        boost::tie(nbIt, nbIt_end) = boost::adjacent_vertices(*ui, fg);
        if (fg[*nbIt].which() == 1) // if the first neighbor is a group it is not a singleton
        {
          //add to set
        }
        else
        {
          //add prot to set
        }
      }
    }
    //sort by score
    //for each:
    // go through peps and remove all incoming connections except the one to this
  }
  */

  void IDBoostGraph::resolveGraphPeptideCentric_(Graph& fg/*, bool resolveTies*/)
  {
    GetPosteriorVisitor gpv{};
    Graph::vertex_iterator ui, ui_end;
    boost::tie(ui,ui_end) = boost::vertices(fg);

    for (; ui != ui_end; ++ui)
    {
      if (fg[*ui].which() == 2)
        // It should suffice to resolve at the pep cluster level
        // if a pep does not belong to a cluster it didnt have multiple parents and
        // therefore does not need to be resolved
      {
        vector<vertex_t> prots;
        queue<vertex_t> start;
        start.push(*ui);
        getUpstreamNodesNonRecursive(start,fg,1,true,prots);
        auto score_compare = [&fg,&gpv](vertex_t& n, vertex_t& m) -> bool
            {return boost::apply_visitor(gpv, fg[n]) < boost::apply_visitor(gpv, fg[m]);};
        auto best_prot = std::max_element(prots.begin(), prots.end(), score_compare); //returns an iterator
        //TODO how to resolve ties

        for (const auto& prot : prots)
        {
          if (prot == *best_prot)
          {
            boost::remove_edge(prot, *ui, fg);
          }
        }
        //TODO remove edges from ID structure as well?
        // if the node is a group, find their members first.
      }
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

     OPENMS_LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;

      #ifdef INFERENCE_MT_DEBUG
     OPENMS_LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
      #endif

      #ifdef INFERENCE_DEBUG
     OPENMS_LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, curr_cc);
     OPENMS_LOG_INFO << "Printed cc " << i << std::endl;
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
       OPENMS_LOG_INFO << "Printing cc " << i << "with intermediate nodes." << std::endl;
        printGraph(LOG_INFO, curr_cc);
       OPENMS_LOG_INFO << "Printed cc " << i << "with intermediate nodes." << std::endl;
        #endif

      }
      else
      {
       OPENMS_LOG_INFO << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  }*/

  //needs run info annotated.
  void IDBoostGraph::clusterIndistProteinsAndPeptidesAndExtendGraph()
  {
    if (nrPrefractionationGroups_ == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Graph not built with run information!");
    }

    /* This should be taken care of. We require that it was built with runinfo here
    if (!pepHitVtx_to_run_.empty()) //graph built with run info
    {
      StringList runs;
      protIDs_.getPrimaryMSRunPath(runs);
      nrReplicates = runs.size();
    }
     */

    pair<int,int> chargeRange = protIDs_.getSearchParameters().getChargeRange();

    if (ccs_.empty()) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "No connected components annotated. Run computeConnectedComponents first!");
    }

    // add_vertex and add_edge not threadsafe
    //#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
    {
      Graph& curr_cc = ccs_[i];

      #ifdef INFERENCE_MT_DEBUG
     OPENMS_LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
      #endif

      #ifdef INFERENCE_DEBUG
     OPENMS_LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;
     OPENMS_LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, curr_cc);
     OPENMS_LOG_INFO << "Printed cc " << i << std::endl;
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
          if (curr_cc[*ui].which() == 0) //protein: same seq peptideHits have to be at a single protein
          {
            SequenceToReplicateChargeVariantHierarchy hierarchy{nrPrefractionationGroups_, chargeRange.first, chargeRange.second};
            Graph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, curr_cc);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              //pepHit, this also makes sure that pepHits already in hierarchy are masked
              if (curr_cc[*adjIt].which() == 6)
              {
                PeptideHit *phitp = boost::get<PeptideHit *>(curr_cc[*adjIt]);
                String seq = phitp->getSequence().toUnmodifiedString();

                //TODO I think it is also best to completely focus on the extended Model here and assume that
                // this information is present. If we allow mixtures of graphical models it gets complex
                // with a lot of if-cases, also/especially during translation to the factor graph.
                Size rep = 0; //In case no replicate info was read.
                if (!pepHitVtx_to_run_.empty()) rep = pepHitVtx_to_run_[*adjIt];
                int chg = phitp->getCharge();

                hierarchy.insert(seq, rep, chg, *adjIt);
              }
            }
            hierarchy.insertToGraph(*ui, curr_cc);
          }
        }

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
              if (curr_cc[*adjIt].which() >= 3) //if there are only two types (pep,prot) this check for pep is actually unnecessary
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
          if (curr_cc[*ui].which() >= 3) //peptide: find peptide clusters
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
       OPENMS_LOG_INFO << "Printing cc " << i << "with intermediate nodes." << std::endl;
        printGraph(LOG_INFO, curr_cc);
       OPENMS_LOG_INFO << "Printed cc " << i << "with intermediate nodes." << std::endl;
        #endif
      }
      else
      {
       OPENMS_LOG_INFO << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  }

  void IDBoostGraph::clusterIndistProteinsAndPeptides()
  {
    if (ccs_.empty()) {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No connected components annotated. Run computeConnectedComponents first!");
    }

    // add_vertex and add_edge not threadsafe
    //#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(ccs_.size()); i += 1)
    {
      Graph& curr_cc = ccs_[i];

      #ifdef INFERENCE_MT_DEBUG
     OPENMS_LOG_INFO << "Processing on thread# " << omp_get_thread_num() << std::endl;
      #endif

      #ifdef INFERENCE_DEBUG
     OPENMS_LOG_INFO << "Processing cc " << i << " with " << boost::num_vertices(curr_cc) << " vertices." << std::endl;
     OPENMS_LOG_INFO << "Printing cc " << i << std::endl;
      printGraph(LOG_INFO, curr_cc);
     OPENMS_LOG_INFO << "Printed cc " << i << std::endl;
      #endif

      // Skip cc without peptide or protein
      //TODO better to do quick bruteforce calculation if the cc is really small
      if (boost::num_edges(curr_cc) >= 1)
      {
        Graph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(curr_cc);

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
              if (curr_cc[*adjIt].which() >= 3) //if there are only two types (pep,prot) this check for pep is actually unnecessary
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
          if (curr_cc[*ui].which() == 6) //peptide: find peptide clusters
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
       OPENMS_LOG_INFO << "Printing cc " << i << "with intermediate nodes." << std::endl;
        printGraph(LOG_INFO, curr_cc);
       OPENMS_LOG_INFO << "Printed cc " << i << "with intermediate nodes." << std::endl;
        #endif
      }
      else
      {
       OPENMS_LOG_INFO << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  }


  //TODO we should probably rename it to splitCC now. Add logging and timing?
  void IDBoostGraph::computeConnectedComponents()
  {
    auto vis = dfs_ccsplit_visitor(ccs_);
    boost::depth_first_search(g, visitor(vis));
   OPENMS_LOG_INFO << "Found " << ccs_.size() << " connected components." << std::endl;
    #ifdef INFERENCE_BENCH
    sizes_and_times_.resize(ccs_.size());
    #endif
    g.clear();
  }

  const IDBoostGraph::Graph& IDBoostGraph::getComponent(Size cc)
  {
    if (cc == 0 && boost::num_vertices(g) != 0)
    {
      return g;
    }
    else
    {
      return ccs_.at(cc);
    }
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

  Size IDBoostGraph::getNrConnectedComponents()
  {
    return ccs_.size();
  }

  ProteinIdentification& IDBoostGraph::getProteinIDs()
  {
    return protIDs_;
  }

  void IDBoostGraph::printGraph(std::ostream& out, const Graph& fg)
  {
    LabelVisitor lv;
    // Also tried to save the labels in a member after build_graph. But could not get the type right for a member that would store them.
    //TODO Is passing "this" to lambda bad? How can I pass private members then?
    auto labels = boost::make_transform_value_property_map([&](const IDPointer &p) { return boost::apply_visitor(lv, p); },
                                                           boost::get(boost::vertex_bundle, fg));
    //boost::print_graph(fg);
    boost::write_graphviz(out, fg, boost::make_label_writer(labels));
  }


  namespace Internal
  {
    /// Hashers for the strong typedefs
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
}

