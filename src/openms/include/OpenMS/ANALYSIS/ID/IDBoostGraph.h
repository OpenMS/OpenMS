// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

// define to get timings for connected components
//#define INFERENCE_BENCH

#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h> //included in BPI
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <unordered_map>
#include <queue>

#include <boost/function.hpp>
#include <boost/blank.hpp>
#include <boost/serialization/strong_typedef.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/properties.hpp>
#include <boost/variant.hpp>
#include <boost/variant/detail/hash_variant.hpp>
#include <boost/variant/static_visitor.hpp>

namespace OpenMS
{
  struct ScoreToTgtDecLabelPairs;

  namespace Internal
  {

  /**
   @brief Creates and maintains a boost graph based on the OpenMS ID datastructures

   For finding connected components and applying functions to them.
   Currently assumes that all PeptideIdentifications are from the ProteinID run that is given.
   Please make sure this is right.
   VERY IMPORTANT NOTE: If you add Visitors here, make sure they do not touch members of the
   underlying ID objects that are responsible for the graph structure. E.g. the (protein/peptide)_hits vectors
   or the lists in ProteinGroups. You can set information like scores or metavalues, though.

   @ingroup Analysis_ID
   */
   //TODO Add OPENMS_DLLAPI everywhere
  class OPENMS_DLLAPI IDBoostGraph
  {

  public:

    // boost has a weird extra semicolon in their strong typedef
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wextra-semi"

    /// placeholder for peptides with the same parent proteins or protein groups
    BOOST_STRONG_TYPEDEF(boost::blank, PeptideCluster);

    /// indistinguishable protein groups (size, nr targets, score)
    struct ProteinGroup
    {
      int size = 0;
      int tgts = 0;
      double score = 0.;
    };

    /// an (currently unmodified) peptide sequence
    BOOST_STRONG_TYPEDEF(String, Peptide);

    /// in which run a PSM was observed
    BOOST_STRONG_TYPEDEF(Size, RunIndex);

    /// in which charge state a PSM was observed
    BOOST_STRONG_TYPEDEF(int, Charge);

    #pragma clang diagnostic pop

    //typedefs
    //TODO rename ProteinGroup type since it collides with the actual OpenMS ProteinGroup
    typedef boost::variant<ProteinHit*, ProteinGroup, PeptideCluster, Peptide, RunIndex, Charge, PeptideHit*> IDPointer;
    typedef boost::variant<const ProteinHit*, const ProteinGroup*, const PeptideCluster*, const Peptide, const RunIndex, const Charge, const PeptideHit*> IDPointerConst;
    //TODO check the impact of different data structures to store nodes/edges
    // Directed graphs would make the internal computations much easier (less in/out edge checking) but boost
    // does not allow computation of "non-strongly" connected components for directed graphs, which is what we would
    // need. We can think about after/while copying to CCs, to insert it into a directed graph!
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer> Graph;
    typedef std::vector<Graph> Graphs;
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer> GraphConst;

    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

    typedef std::set<IDBoostGraph::vertex_t> ProteinNodeSet;
    typedef std::set<IDBoostGraph::vertex_t> PeptideNodeSet;


    /// A boost dfs visitor that copies connected components into a vector of graphs
    class dfs_ccsplit_visitor:
        public boost::default_dfs_visitor
    {
    public:
      dfs_ccsplit_visitor(Graphs& vgs)
          :  gs(vgs), curr_v(0), next_v(0), m()
      {}

      template < typename Vertex, typename Graph >
      void start_vertex(Vertex u, const Graph & tg)
      {
        gs.emplace_back();
        next_v = boost::add_vertex(tg[u], gs.back());
        m[u] = next_v;
      }

      template < typename Vertex, typename Graph >
      void discover_vertex(Vertex /*u*/, const Graph & /*tg*/)
      {
        curr_v = next_v;
      }

      template < typename Edge, typename Graph >
      void examine_edge(Edge e, const Graph & tg)
      {
        if (m.find(e.m_target) == m.end())
        {
          next_v = boost::add_vertex(tg[e.m_target], gs.back());
          m[e.m_target] = next_v;
        }
        else
        {
          next_v = m[e.m_target];
        }

        boost::add_edge(m[e.m_source], next_v, gs.back());
      }

      Graphs& gs;
      vertex_t curr_v, next_v;
      /// A mapping from old node id to new node id to not duplicate existing ones in the new graph
      std::map<vertex_t, vertex_t> m;
    };

    ///@brief Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type creates a label
    /// e.g. for printing to dot format
    class LabelVisitor:
        public boost::static_visitor<OpenMS::String>
    {
    public:

      OpenMS::String operator()(const PeptideHit* pep) const
      {
        return pep->getSequence().toString() + "_" + pep->getCharge();
      }

      OpenMS::String operator()(const ProteinHit* prot) const
      {
        return prot->getAccession();
      }

      OpenMS::String operator()(const ProteinGroup& /*protgrp*/) const
      {
        return "PG";
      }

      OpenMS::String operator()(const PeptideCluster& /*pc*/) const
      {
        return "PepClust";
      }

      OpenMS::String operator()(const Peptide& peptide) const
      {
        return peptide;
      }

      OpenMS::String operator()(const RunIndex& ri) const
      {
        return "rep" + String(ri);
      }

      OpenMS::String operator()(const Charge& chg) const
      {
        return "chg" + String(chg);
      }

    };

    /// @brief Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type prints the address.
    /// For debugging purposes only
    template<class CharT>
    class PrintAddressVisitor:
        public boost::static_visitor<>
    {
    public:

      explicit PrintAddressVisitor(std::basic_ostream<CharT> stream):
          stream_(stream)
      {}

      void operator()(PeptideHit* pep) const
      {
        stream_ << pep->getSequence().toUnmodifiedString() << ": " << pep << std::endl;
      }

      void operator()(ProteinHit* prot) const
      {
        stream_ << prot->getAccession() << ": " << prot << std::endl;
      }

      void operator()(const ProteinGroup& /*protgrp*/) const
      {
        stream_ << "PG" << std::endl;
      }

      void operator()(const PeptideCluster& /*pc*/) const
      {
        stream_ << "PepClust" << std::endl;
      }

      void operator()(const Peptide& peptide) const
      {
        stream_ << peptide << std::endl;
      }

      void operator()(const RunIndex& ri) const
      {
        stream_ << "rep" << ri << std::endl;
      }

      void operator()(const Charge& chg) const
      {
        stream_ << "chg" << chg << std::endl;
      }

      std::basic_ostream<CharT> stream_;
    };

    /// @brief Visits nodes in the boost graph (either ptrs to an ID Object or some lightweight surrogates)
    /// and depending on their type sets the posterior
    /// Don't forget to set higherScoreBetter and score names in the parent ID objects.
    class SetPosteriorVisitor:
        public boost::static_visitor<>
    {
    public:

      void operator()(PeptideHit* pep, double posterior) const
      {
        pep->setScore(posterior);
      }

      void operator()(ProteinHit* prot, double posterior) const
      {
        prot->setScore(posterior);
      }

      void operator()(ProteinGroup& pg, double posterior) const
      {
        pg.score = posterior;
      }

      // Everything else, do nothing for now
      template <class T>
      void operator()(T& /*any node type*/, double /*posterior*/) const
      {
        // do nothing
      }

    };

    /// @brief Visits nodes in the boost graph (either ptrs to an ID Object or some lightweight surrogates)
    /// and depending on their type gets the score (usually the posterior)
    class GetPosteriorVisitor:
        public boost::static_visitor<double>
    {
    public:

      double operator()(PeptideHit* pep) const
      {
        return pep->getScore();
      }

      double operator()(ProteinHit* prot) const
      {
        return prot->getScore();
      }

      double operator()(ProteinGroup& pg) const
      {
        return pg.score;
      }

      // Everything else, do nothing for now
      template <class T>
      double operator()(T& /*any node type*/) const
      {
        return -1.0;
      }

    };

    /// @brief Visits nodes in the boost graph (either ptrs to an ID Object or some lightweight surrogates)
    /// and depending on their type gets the score (usually the posterior) plus if it is a decoy or a target.
    /// If not known or not defined, returns (-1.0, false)
    class GetScoreTgTVisitor:
    public boost::static_visitor<std::pair<double,bool>>
        {
        public:

          std::pair<double,bool> operator()(PeptideHit* pep) const
          {
            return {pep->getScore(), pep->getMetaValue("target_decoy").toString()[0] == 't'};
          }

          std::pair<double,bool> operator()(ProteinHit* prot) const
          {
            return {prot->getScore(), prot->getMetaValue("target_decoy").toString()[0] == 't'};
          }

          std::pair<double,bool> operator()(ProteinGroup& pg) const
          {
            return {pg.score, pg.tgts > 0};
          }

          // Everything else, do nothing for now
          template <class T>
          std::pair<double,bool> operator()(T& /*any node type*/) const
          {
            return {-1.0, false};
          }
    };

    /// Constructors
    IDBoostGraph(ProteinIdentification& proteins,
                std::vector<PeptideIdentification>& idedSpectra,
                Size use_top_psms,
                bool use_run_info,
                bool best_psms_annotated,
                const std::optional<const ExperimentalDesign>& ed = std::optional<const ExperimentalDesign>());

    IDBoostGraph(ProteinIdentification& proteins,
                ConsensusMap& cmap,
                Size use_top_psms,
                bool use_run_info,
                bool use_unassigned_ids,
                bool best_psms_annotated,
                const std::optional<const ExperimentalDesign>& ed = std::optional<const ExperimentalDesign>());


    //TODO think about templating to avoid wrapping to std::function
    // although we usually do long-running tasks per CC such that the extra virtual call does not matter much
    // Instead we gain type erasure.
    /// Do sth on connected components (your functor object has to inherit from std::function or be a lambda)
    void applyFunctorOnCCs(const std::function<unsigned long(Graph&, unsigned int)>& functor);
    /// Do sth on connected components single threaded (your functor object has to inherit from std::function or be a lambda)
    void applyFunctorOnCCsST(const std::function<void(Graph&)>& functor);

    /// Add intermediate nodes to the graph that represent indist. protein groups and peptides with the same parents
    /// this will save computation time and oscillations later on.
    void clusterIndistProteinsAndPeptides();

    //TODO create a new class for an extended Graph and try to reuse as much as possible
    // use inheritance or templates
    /// (under development) As above but adds charge, replicate and sequence layer of nodes (untested)
    /// @todo needs to be finished, updated with latest additions (i.e. check clusterIndistProteinsAndPeptides), and tested
    void clusterIndistProteinsAndPeptidesAndExtendGraph();

    /// Annotate indistinguishable proteins by adding the groups to the underlying
    /// ProteinIdentification::ProteinGroups object.
    /// This has no effect on the graph itself.
    /// @pre Graph must contain ProteinGroup nodes (e.g. with clusterIndistProteinsAndPeptides).
    /// Otherwise it does nothing and you should use calculateAndAnnotateIndistProteins instead.
    /// @param addSingletons if you want to annotate groups with just one protein entry
    void annotateIndistProteins(bool addSingletons = true);

    /// Annotate indistinguishable proteins by adding the groups to the underlying
    /// ProteinIdentification::ProteinGroups object. This has no effect on the graph itself.
    /// @param addSingletons if you want to annotate groups with just one protein entry
    void calculateAndAnnotateIndistProteins(bool addSingletons = true);

    /// Splits the initialized graph into connected components and clears it.
    void computeConnectedComponents();

    /// @todo untested
    /// Removes all edges from a peptide (and its PSMs) to its parent protein groups (and its proteins)
    /// except for the best protein group.
    /// @pre Graph must contain PeptideCluster nodes (e.g. with clusterIndistProteinsAndPeptides).
    /// @param removeAssociationsInData Also removes the corresponding PeptideEvidences in the underlying
    ///     ID data structure. Only deactivate if you know what you are doing.
    void resolveGraphPeptideCentric(bool removeAssociationsInData = true);



    /// Zero means the graph was not split yet
    Size getNrConnectedComponents();

    /// @brief Returns a specific connected component of the graph as a graph itself
    /// @param cc the index of the component
    /// @return the component as graph
    const Graph& getComponent(Size cc);

    /// @brief Returns the underlying protein identifications for viewing
    /// @return const ref to the protein ID run in this graph (can only be one)
    const ProteinIdentification& getProteinIDs();

    //TODO docu
    //void buildExtendedGraph(bool use_all_psms, std::pair<int,int> chargeRange, unsigned int nrReplicates);

    /// @brief Prints a graph (component or if not split, the full graph) in graphviz (i.e. dot) format
    /// @param out an ostream to print to
    /// @param fg the graph to print
    static void printGraph(std::ostream& out, const Graph& fg);

    /// @brief Searches for all upstream nodes from a (set of) start nodes that are lower
    ///    or equal than a given level. The ordering is the same as in the IDPointer variant typedef.
    /// @param q a queue of start nodes
    /// @param graph the graph to look in (q has to be part of it)
    /// @param lvl the level to start reporting from
    /// @param stop_at_first do you want to stop at the first node <= lvl or also report its
    ///    upstream "predecessors"
    /// @param result vector of reported nodes
    void getUpstreamNodesNonRecursive(std::queue<vertex_t>& q, const Graph& graph, int lvl,
                                      bool stop_at_first, std::vector<vertex_t>& result);

    /// @brief Searches for all downstream nodes from a (set of) start nodes that are higher
    ///    or equal than a given level. The ordering is the same as in the IDPointer variant typedef.
    /// @param q a queue of start nodes
    /// @param graph the graph to look in (q has to be part of it)
    /// @param lvl the level to start reporting from
    /// @param stop_at_first do you want to stop at the first node >= lvl or also report its
    ///    upstream "predecessors"
    /// @param result vector of reported nodes
    void getDownstreamNodesNonRecursive(std::queue<vertex_t>& q, const Graph& graph, int lvl,
                                        bool stop_at_first, std::vector<vertex_t>& result);

    /// Gets the scores from the proteins included in the graph.
    /// The difference to querying the underlying ProteinIdentification structure is that not all
    /// proteins might be included in the graph due to using only the best psm per peptide
    void getProteinScores_(ScoreToTgtDecLabelPairs& scores_and_tgt);
    /// Gets the scores and target decoy fraction from groups and score + binary values for singleton
    /// proteins. This function is usually used to create input for FDR calculations
    void getProteinGroupScoresAndTgtFraction(ScoreToTgtDecLabelPairs& scores_and_tgt_fraction);
    void getProteinGroupScoresAndHitchhikingTgtFraction(ScoreToTgtDecLabelPairs& scores_and_tgt_fraction);

  private:

    ProteinIdentification& protIDs_;

    struct SequenceToReplicateChargeVariantHierarchy;


    //TODO introduce class hierarchy:
    /*
     * IDGraph<UnderlyingIDStruc>
     *
     * - BasicGraph<>
     * - ExtendedGraphClustered<>
     * - ExtendedGraphClusteredWithRunInfo<>
     *
     * in theory extending a basic one is desirable to create the extended one. But it means we have to
     * copy/move the graph (node by node) because the nodes are of a broader boost::variant type. So we probably have to
     * duplicate code and offer a from-scratch step-wise building for the extended graph, too.
     * Note that there could be several levels of extension in the future. For now I keep everything in one
     * class by having potential storage for the broadest extended type. Differences in the underlying ID structure
     * e.g. ConsensusMap or PeptideIDs from idXML currently only have an effect during building, so I just overload
     * the constructors. In theory it would be nice to generalize on that, too, especially when we adapt to the new
     * ID data structure.
     */


    /* ----------------  Either of them is used, preferably second  --------------- */
    /// the initial boost Graph (will be cleared when split into CCs)
    Graph g;

    /// the Graph split into connected components
    Graphs ccs_;
    /* ---------------------------------------------------------------------------- */

    #ifdef INFERENCE_BENCH
    /// nrNodes, nrEdges, nrMessages and times of last functor execution per connected component
    std::vector<std::tuple<vertex_t, vertex_t, unsigned long, double>> sizes_and_times_{1};
    #endif


    /* ----  Only used when run information was available --------- */

    //TODO think about preallocating it, but the number of peptide hits is not easily computed
    // since they are inside the pepIDs
    //TODO would multiple sets be better?

    /// if a graph is built with run information, this will store the run, each peptide hit
    /// vertex belongs to. Important for extending the graph.
    std::unordered_map<vertex_t, Size> pepHitVtx_to_run_;

    /// this basically stores the number of different values in the pepHitVtx_to_run
    /// a Prefractionation group (previously called run) is a unique combination
    /// of all non-fractionation related entries in the exp. design
    /// i.e. one (sub-)experiment before fractionation
    Size nrPrefractionationGroups_ = 0;

    /* ----------------------------------------------------------- */


    /// helper function to add a vertex if it is not present yet, otherwise return the present one
    /// needs a temporary filled vertex_map that is modifiable
    vertex_t addVertexWithLookup_(const IDPointer& ptr, std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map);
    //vertex_t addVertexWithLookup_(IDPointerConst& ptr, std::unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>>& vertex_map);


    /// internal function to annotate the underlying ID structures based on the given Graph
    void annotateIndistProteins_(const Graph& fg, bool addSingletons);
    void calculateAndAnnotateIndistProteins_(const Graph& fg, bool addSingletons);

    /// Initialize and store the graph
    /// IMPORTANT: Once the graph is built, editing members like (protein/peptide)_hits_ will invalidate it!
    /// @param proteins ProteinIdentification object storing IDs and groups
    /// @param idedSpectra vector of ProteinIdentifications with links to the proteins and PSMs in its PeptideHits
    /// @param use_top_psms Nr of top PSMs used per spectrum (<= 0 means all)
    /// @param best_psms_annotated Are the PSMs annotated with the "best_per_peptide" meta value. Otherwise all are
    ///  taken into account.
    /// @todo we could include building the graph in important "main" functions like inferPosteriors
    /// to make the methods safer, but it is also nice to be able to reuse the graph
    void buildGraph_(ProteinIdentification& proteins,
                    std::vector<PeptideIdentification>& idedSpectra,
                    Size use_top_psms,
                    bool best_psms_annotated = false);

    void buildGraph_(ProteinIdentification& proteins,
                    ConsensusMap& cmap,
                    Size use_top_psms,
                    bool use_unassigned_ids,
                    bool best_psms_annotated = false);

    /// Used during building
    void addPeptideIDWithAssociatedProteins_(
        PeptideIdentification& spectrum,
        std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map,
        const std::unordered_map<std::string, ProteinHit*>& accession_map,
        Size use_top_psms,
        bool best_psms_annotated);

    void addPeptideAndAssociatedProteinsWithRunInfo_(
        PeptideIdentification& spectrum,
        std::unordered_map<unsigned, unsigned>& indexToPrefractionationGroup,
        std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map,
        std::unordered_map<std::string, ProteinHit*>& accession_map,
        Size use_top_psms
        );

    /// Initialize and store the graph. Also stores run information to later group
    /// peptides more efficiently.
    /// IMPORTANT: Once the graph is built, editing members like (protein/peptide)_hits_ will invalidate it!
    /// 
    /// @p use_top_psms is the number of top PSMs used per spectrum (<= 0 means all)
    /// @todo we could include building the graph in important "main" functions like inferPosteriors
    /// to make the methods safer, but it is also nice to be able to reuse the graph
    void buildGraphWithRunInfo_(ProteinIdentification& proteins,
                               ConsensusMap& cmap,
                               Size use_top_psms,
                               bool use_unassigned_ids,
                               const ExperimentalDesign& ed);

    void buildGraphWithRunInfo_(ProteinIdentification& proteins,
                               std::vector<PeptideIdentification>& idedSpectra,
                               Size use_top_psms,
                               const ExperimentalDesign& ed);


    /// see equivalent public method
    void resolveGraphPeptideCentric_(Graph& fg, bool removeAssociationsInData);

    template<class NodeType>
    void getDownstreamNodes(const vertex_t& start, const Graph& graph, std::vector<NodeType>& result)
    {
      Graph::adjacency_iterator adjIt, adjIt_end;
      boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(start, graph);
      for (;adjIt != adjIt_end; ++adjIt)
      {
        if (graph[*adjIt].type() == typeid(NodeType))
        {
          result.emplace_back(boost::get<NodeType>(graph[*adjIt]));
        }
        else if (graph[*adjIt].which() > graph[start].which())
        {
          getDownstreamNodes(*adjIt, graph, result);
        }
      }
    }

    template<class NodeType>
    void getUpstreamNodes(const vertex_t& start, const Graph graph, std::vector<NodeType>& result)
    {
      Graph::adjacency_iterator adjIt, adjIt_end;
      boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(start, graph);
      for (;adjIt != adjIt_end; ++adjIt)
      {
        if (graph[*adjIt].type() == typeid(NodeType))
        {
          result.emplace_back(boost::get<NodeType>(graph[*adjIt]));
        }
        else if (graph[*adjIt].which() < graph[start].which())
        {
          getUpstreamNodes(*adjIt, graph, result);
        }
      }
    }
  };

    bool operator==(const IDBoostGraph::ProteinGroup& lhs, const IDBoostGraph::ProteinGroup& rhs);
  } //namespace Internal
} //namespace OpenMS

