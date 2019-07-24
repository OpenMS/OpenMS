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

#pragma once

// define to get timings for connected components
//#define INFERENCE_BENCH

#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h> //included in BPI
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <unordered_map>
#include <queue>

#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/properties.hpp>
#include <boost/variant.hpp>
#include <boost/variant/detail/hash_variant.hpp>
#include <boost/variant/static_visitor.hpp>

namespace OpenMS {
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
    BOOST_STRONG_TYPEDEF(boost::blank, PeptideCluster)

    /// indistinguishable protein groups
    BOOST_STRONG_TYPEDEF(double, ProteinGroup)

    /// an (currently unmodified) peptide sequence
    BOOST_STRONG_TYPEDEF(String, Peptide)

    /// in which run a PSM was observed
    BOOST_STRONG_TYPEDEF(Size, RunIndex)

    /// in which charge state a PSM was observed
    BOOST_STRONG_TYPEDEF(int, Charge)

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

    //TODO group visitors by templates
    /// Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type creates a label
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
        return String("PG");
      }

      OpenMS::String operator()(const PeptideCluster& /*pc*/) const
      {
        return String("PepClust");
      }

      OpenMS::String operator()(const Peptide& peptide) const
      {
        return peptide;
      }

      OpenMS::String operator()(const RunIndex& ri) const
      {
        return String("rep" + String(ri));
      }

      OpenMS::String operator()(const Charge& chg) const
      {
        return String("chg" + String(chg));
      }

    };

    /// Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type prints the address.
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

    /// Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type sets the posterior
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
        pg = posterior;
      }

      // Everything else, do nothing for now
      template <class T>
      void operator()(T& /*any node type*/, double /*posterior*/) const
      {
        // do nothing
      }

    };

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
        return pg;
      }

      // Everything else, do nothing for now
      template <class T>
      double operator()(T& /*any node type*/) const
      {
        return -1.0;
      }

    };

    /// Constructors
    IDBoostGraph(ProteinIdentification& proteins,
                std::vector<PeptideIdentification>& idedSpectra,
                Size use_top_psms,
                bool use_run_info,
                const boost::optional<const ExperimentalDesign>& ed = boost::optional<const ExperimentalDesign>());

    IDBoostGraph(ProteinIdentification& proteins,
                ConsensusMap& cmap,
                Size use_top_psms,
                bool use_run_info,
                bool use_unassigned_ids,
                const boost::optional<const ExperimentalDesign>& ed = boost::optional<const ExperimentalDesign>());


    /// Do sth on connected components (your functor object has to inherit from std::function or be a lambda)
    void applyFunctorOnCCs(const std::function<unsigned long(Graph&)>& functor);
    /// Do sth on connected components single threaded (your functor object has to inherit from std::function or be a lambda)
    void applyFunctorOnCCsST(const std::function<void(Graph&)>& functor);

    /// Add intermediate nodes to the graph that represent indist. protein groups and peptides with the same parents
    /// this will save computation time and oscillations later on.
    void clusterIndistProteinsAndPeptides();

    //TODO create a new class for an extended Graph and try to reuse as much as possible
    // use inheritance or templates
    /// (under development) As above but adds charge, replicate and sequence layer of nodes (untested)
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



    /// Zero means the graph was not split yet
    Size getNrConnectedComponents();

    const Graph& getComponent(Size cc);

    ProteinIdentification& getProteinIDs();

    //TODO docu
    //void buildExtendedGraph(bool use_all_psms, std::pair<int,int> chargeRange, unsigned int nrReplicates);

    static void printGraph(std::ostream& out, const Graph& fg);

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
    /// nrnodes, nredges, nrmessages and times of last functor execution per connected component
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
    vertex_t addVertexWithLookup_(IDPointer& ptr, std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map);
    //vertex_t addVertexWithLookup_(IDPointerConst& ptr, std::unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>>& vertex_map);


    /// internal function to annotate the underlying ID structures based on the given Graph
    void annotateIndistProteins_(const Graph& fg, bool addSingletons);
    void calculateAndAnnotateIndistProteins_(const Graph& fg, bool addSingletons);

    /// Initialize and store the graph
    /// IMPORTANT: Once the graph is built, editing members like (protein/peptide)_hits_ will invalidate it!
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param idedSpectra vector of ProteinIdentifications with links to the proteins and PSMs in its PeptideHits
    /// @param use_top_psms Nr of top PSMs used per spectrum (<= 0 means all)
    /// @todo we could include building the graph in important "main" functions like inferPosteriors
    /// to make the methods safer, but it is also nice to be able to reuse the graph
    void buildGraph_(ProteinIdentification& proteins,
                    std::vector<PeptideIdentification>& idedSpectra,
                    Size use_top_psms);

    void buildGraph_(ProteinIdentification& proteins,
                    ConsensusMap& cmap,
                    Size use_top_psms,
                    bool use_unassigned_ids);

    /// Used during building
    void addPeptideIDWithAssociatedProteins_(
        PeptideIdentification& spectrum,
        std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map,
        const std::unordered_map<std::string, ProteinHit*>& accession_map,
        Size use_top_psms);

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
    /// @param use_top_psms Nr of top PSMs used per spectrum (<= 0 means all)
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


    void getUpstreamNodesNonRecursive(std::queue<vertex_t>& q, Graph graph, int lvl,
        bool stop_at_first, std::vector<vertex_t>& result);

    void resolveGraphPeptideCentric_(Graph& fg);

    template<class NodeType>
    void getDownstreamNodes(vertex_t start, Graph graph, std::vector<NodeType>& result)
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
    void getUpstreamNodes(vertex_t start, Graph graph, std::vector<NodeType>& result)
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

} } //namespace OpenMS

