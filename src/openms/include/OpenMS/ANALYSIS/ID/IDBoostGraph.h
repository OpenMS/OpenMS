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

#ifndef OPENMS_ANALYSIS_ID_IDBOOSTGRAPH_H
#define OPENMS_ANALYSIS_ID_IDBOOSTGRAPH_H

// define to get timings for connected components
#define INFERENCE_BENCH

#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h> //included in BPI
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <unordered_map>

#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/properties.hpp>
#include <boost/variant.hpp>
#include <boost/variant/detail/hash_variant.hpp>
#include <boost/variant/static_visitor.hpp>

namespace OpenMS
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

    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wextra-semi"

    BOOST_STRONG_TYPEDEF(boost::blank, PeptideCluster)

    BOOST_STRONG_TYPEDEF(double, ProteinGroup)

    BOOST_STRONG_TYPEDEF(String, Peptide)

    BOOST_STRONG_TYPEDEF(Size, RunIndex)

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
    //typedef boost::subgraph<boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer>> SubGraph;
    typedef std::vector<Graph> Graphs;
    //typedef std::vector<SubGraph> SubGraphs;
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer> GraphConst;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
    typedef boost::filtered_graph<Graph, boost::function<bool(edge_t)>, boost::function<bool(vertex_t)> > FilteredGraph;
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
        #ifdef INFERENCE_DEBUG
        std::cout << "set score " << posterior << " for " << prot->getAccession() << std::endl;
        #endif
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

    /// Constructors
    IDBoostGraph(ProteinIdentification &proteins, std::vector<PeptideIdentification>& idedSpectra);
    IDBoostGraph(
        ProteinIdentification &proteins,
        std::vector<PeptideIdentification>& idedSpectra,
        const ExperimentalDesign& ed);

    /// Do sth on connected components (your functor object has to inherit from std::function)
    void applyFunctorOnCCs(std::function<unsigned long(Graph&)> functor);
    void applyFunctorOnCCsST(std::function<void(Graph&)> functor);

    /// Add intermediate nodes to the graph that represent indist. protein groups and peptides with the same parents
    /// this will save computation time and oscillations later on.
    void clusterIndistProteinsAndPeptides();
    /// As above but adds charge, replicate and sequence layer of nodes (untested)
    void clusterIndistProteinsAndPeptidesAndExtendGraph();

    /// Annotate indistinguishable proteins by adding the groups to the underlying
    /// ProteinIdentification::ProteinGroups object. This has no effect on the graph itself.
    /// @param addSingletons if you want to annotate groups with just one protein entry
    void annotateIndistProteins(bool addSingletons = true) const;

    /// Splits the initialized graph into connected components and clears it.
    void computeConnectedComponents();

    /// Initialize and store the graph
    /// IMPORTANT: Once the graph is built, editing members like (protein/peptide)_hits_ will invalidate it!
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param idedSpectra vector of ProteinIdentifications with links to the proteins and PSMs in its PeptideHits
    /// @param use_all_psms If all or just the FIRST psm should be used
    void buildGraph(Size use_top_psms);
    void buildGraphWithRunInfo(Size use_top_psms, bool readstore_run_info = true);

    Size getNrConnectedComponents();

    //TODO docu
    //void buildExtendedGraph(bool use_all_psms, std::pair<int,int> chargeRange, unsigned int nrReplicates);

    static void printFilteredGraph(std::ostream& out, const FilteredGraph& fg);
    static void printGraph(std::ostream& out, const Graph& fg);

  private:

    struct SequenceToReplicateChargeVariantHierarchy;

    /// the initial boost Graph
    Graph g;

    /// the Graph split into connected components
    Graphs ccs_;

    #ifdef INFERENCE_BENCH
    /// nrnodes, nredges, nrmessages and times of last functor execution per connected component
    std::vector<std::tuple<vertex_t, vertex_t, unsigned long, double>> sizes_and_times_{1};
    #endif

    /// underlying protein identification object
    //TODO for consensusXML this probably needs to become a vector.
    ProteinIdentification& proteins_;
    /// underlying peptide identifications
    std::vector<PeptideIdentification>& idedSpectra_;
    /// underlying experimental design, if not given it will be default constructed
    //TODO think about using a pointer here instead of copying. But it's usually not too big
    const ExperimentalDesign exp_design_;

    /// if a graph is built with run information, this will store the run, each peptide hit
    /// vertex belongs to. Important for extending the graph.
    //TODO think about preallocating it, but the number of peptide hits is not easily computed
    //since they are inside the pepIDs
    //TODO would multiple sets be better?
    std::unordered_map<vertex_t, Size> pepHitVtx_to_run_;

    /// a visitor that creates labels based on the node type (e.g. for printing)
    //LabelVisitor lv_; //currently created locally

    /// helper function to add a vertex if it is not present yet, otherwise return the present one
    /// needs a temporary filled vertex_map that is modifiable
    vertex_t addVertexWithLookup_(IDPointer& ptr, std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map);
    //vertex_t addVertexWithLookup_(IDPointerConst& ptr, std::unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>>& vertex_map);


    /// internal function to annotate the underlying ID structures based on the given Graph
    void annotateIndistProteins_(const Graph& fg, bool addSingletons) const;

  };

} //namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDBOOSTGRAPH_H
