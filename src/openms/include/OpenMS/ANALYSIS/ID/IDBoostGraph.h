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

#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h> //included in BPI
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <unordered_map>

#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/subgraph.hpp>
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
  class IDBoostGraph
  {

  public:

    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wextra-semi"

    BOOST_STRONG_TYPEDEF(char, PeptideCluster)
    std::size_t hash_value(PeptideCluster const& x)
    {
      boost::hash<char> hasher;
      return hasher(static_cast<char>(x));
    }

    // TODO make a double out of it and store the posterior
    BOOST_STRONG_TYPEDEF(char, ProteinGroup)
    std::size_t hash_value(ProteinGroup const& x)
    {
      boost::hash<char> hasher;
      return hasher(static_cast<char>(x));
    }

    BOOST_STRONG_TYPEDEF(String, Peptide)
    std::size_t hash_value(Peptide const& x)
    {
      boost::hash<std::string> hasher;
      return hasher(static_cast<std::string>(x));
    }

    BOOST_STRONG_TYPEDEF(Size, RunIndex)
    std::size_t hash_value(RunIndex const& x)
    {
      boost::hash<Size> hasher;
      return hasher(static_cast<Size>(x));
    }

    BOOST_STRONG_TYPEDEF(int, Charge)
    std::size_t hash_value(Charge const& x)
    {
      boost::hash<int> hasher;
      return hasher(static_cast<int>(x));
    }

    #pragma clang diagnostic pop

    //typedefs
    //typedef ProteinIdentification::ProteinGroup ProteinGroup;

    typedef boost::variant<ProteinHit*, ProteinGroup*, PeptideCluster*, Peptide, RunIndex, Charge, PeptideHit*> IDPointer;
    typedef boost::variant<const ProteinHit*, const ProteinGroup*, const PeptideCluster*, const Peptide, const RunIndex, const Charge, const PeptideHit*> IDPointerConst;
    //TODO check the impact of different data structures to store nodes/edges
    // Directed graphs would make the internal computations much easier (less in/out edge checking) but boost
    // does not allow computation of "non-strongly" connected components for directed graphs, which is what we would
    // need
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer> Graph;
    typedef boost::subgraph<boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer>> SubGraph;
    typedef std::vector<Graph> Graphs;
    typedef std::vector<SubGraph> SubGraphs;
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer> GraphConst;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
    typedef boost::filtered_graph<Graph, boost::function<bool(edge_t)>, boost::function<bool(vertex_t)> > FilteredGraph;
    typedef std::set<IDBoostGraph::vertex_t> ProteinNodeSet;
    typedef std::set<IDBoostGraph::vertex_t> PeptideNodeSet;

    /// Constructor
    IDBoostGraph(ProteinIdentification &proteins, std::vector<PeptideIdentification>& idedSpectra);

    /// Do sth on connected components (your functor object has to inherit from std::function)
    void applyFunctorOnCCs(std::function<void(Graph&)> functor);

    /// Add intermediate nodes to the graph that represent indist. protein groups and peptides with the same parents
    /// this will save computation time and oscillations later on.
    void clusterIndistProteinsAndPeptides();

    /// Annotate indistinguishable proteins by adding the groups to the underlying
    /// ProteinIdentification::ProteinGroups object. This has no effect on the graph itself.
    /// @param addSingletons if you want to annotate groups with just one protein entry
    void annotateIndistProteins(bool addSingletons = true) const;

    //TODO move to cpp file
    struct SequenceToReplicateChargeVariantHierarchy
    {
      //TODO only add the intermediate nodes if there are more than one "splits"
      SequenceToReplicateChargeVariantHierarchy(Size nrReplicates, int minCharge, int maxCharge):
          seq_to_vecs_{},
          minCharge_(minCharge),
          nrCharges_(Size(maxCharge - minCharge) + 1u),
          nrReplicates_(nrReplicates)
      {}

      void insert(String& seq, Size replicate, Size charge, vertex_t pepVtx)
      {
        auto seq_it = seq_to_vecs_.emplace(std::move(seq), std::vector<std::vector<std::set<vertex_t>>>{nrReplicates_, {nrCharges_, {}}});
        seq_it.first->second[replicate][charge - minCharge_].insert(pepVtx);
      }

      void insertToGraph(vertex_t rootProteinVtx, Graph& graph)
      {
        for (const auto& seqContainer : seq_to_vecs_)
        {
          vertex_t pep = boost::add_vertex(Peptide{seqContainer.first}, graph);
          boost::add_edge(rootProteinVtx, pep, graph);
          for (Size s = 0; s < seqContainer.second.size(); ++s)
          {
            //TODO Gather prots for this sequence by getting all proteins attached to one of the decendents
            //std::vector<vertex_t> prots_for_pepseq;
            // for adjacency iterator: if protein : push_back index
            vertex_t ri = boost::add_vertex(RunIndex{s},graph);
            boost::add_edge(pep, ri, graph);
            for (Size t = 0; t < seqContainer.second[s].size(); ++t)
            {
              vertex_t cs = boost::add_vertex(Charge{minCharge_ + int(t)}, graph);
              boost::add_edge(ri, cs, graph);
              for (const auto& pepVtx : seqContainer.second[s][t])
              {
                //TODO ACTUALLY REWIRE EDGES TO ALL PROTEIN PARENTS
                //for parent : prots_for_pepseq, do the things below
                boost::add_edge(cs, pepVtx, graph);
                boost::remove_edge(rootProteinVtx, pepVtx, graph);
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

      OpenMS::String operator()(const ProteinGroup* /*protgrp*/) const
      {
        return String("PG");
      }

      OpenMS::String operator()(const PeptideCluster* /*pc*/) const
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

      void operator()(const ProteinGroup* /*protgrp*/) const
      {
        stream_ << "PG" << std::endl;
      }

      void operator()(const PeptideCluster* /*pc*/) const
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

      void operator()(ProteinGroup* /*protgrp*/, double /*posterior*/) const
      {
        // do nothing
        // TODO we could store posterior here
        // protgrp->probability = posterior;
      }

      void operator()(const PeptideCluster* /*pc*/, double /*posterior*/) const
      {
        // do nothing
      }

      void operator()(const Peptide&, double) const
      {
        // do nothing
      }

      void operator()(const RunIndex&, double) const
      {
        // do nothing
      }

      void operator()(const Charge&, double) const
      {
        // do nothing
      }

    };

    /// Splits the initialized graph into connected components and clears it.
    void computeConnectedComponents();

    /// Initialize and store the graph
    /// IMPORTANT: Once the graph is built, editing members like (protein/peptide)_hits_ will invalidate it!
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param idedSpectra vector of ProteinIdentifications with links to the proteins and PSMs in its PeptideHits
    /// @param use_all_psms If all or just the FIRST psm should be used
    void buildGraph(bool use_all_psms);

    //TODO docu
    //void buildExtendedGraph(bool use_all_psms, std::pair<int,int> chargeRange, unsigned int nrReplicates);

  private:

    /// the initial boost Graph
    Graph g;

    /// the Graph split into connected components
    Graphs ccs_;

    /// static objects created from a hard typedef to differentiate between different types of nodes.
    /// they will be reused throughout the graph when necessary
    static PeptideCluster staticPC;
    static ProteinGroup staticPG;

    //GraphConst gconst;

    /// underlying protein identification object
    //TODO for consensusXML this probably needs to become a vector.
    ProteinIdentification& proteins_;
    /// underlying peptide identifications
    std::vector<PeptideIdentification>& idedSpectra_;

    /// if a graph is built with run information, this will store the run, each peptide hit
    /// vertex belongs to. Important for extending the graph.
    //TODO think about preallocating it, but the number of peptide hits is not easily computed
    //since they are inside the pepIDs
    //TODO would multiple sets be better?
    std::unordered_map<vertex_t, Size> pepHitVtx_to_run_;

    /// a visitor that creates labels based on the node type (e.g. for printing)
    LabelVisitor lv_;

    /// helper function to add a vertex if it is not present yet, otherwise return the present one
    /// needs a temporary filled vertex_map that is modifiable
    vertex_t addVertexWithLookup_(IDPointer& ptr, std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map);
    //vertex_t addVertexWithLookup_(IDPointerConst& ptr, std::unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>>& vertex_map);


    /// internal function to annotate the underlying ID structures based on the given Graph
    void annotateIndistProteins_(const Graph& fg, bool addSingletons) const;

    void printFilteredGraph(std::ostream& out, const FilteredGraph& fg) const;
    void printGraph(std::ostream& out, const Graph& fg) const;

  };

} //namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDBOOSTGRAPH_H
