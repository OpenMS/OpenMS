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

#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
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
   VERY IMPORTANT NOTE: If you add Visitors here, make sure they do not touch members of the
   underlying ID objects that are responsible for the graph structure. E.g. the (protein/peptide)_hits vectors
   or the lists in ProteinGroups. You can set information like scores or metavalues, though.

   @ingroup Analysis_ID
   */
  class IDBoostGraph
  {

  public:

    BOOST_STRONG_TYPEDEF(char, PeptideCluster)
    BOOST_STRONG_TYPEDEF(char, ProteinGroup)

    //typedefs
    //typedef ProteinIdentification::ProteinGroup ProteinGroup;

    typedef boost::variant<ProteinHit*, ProteinGroup*, PeptideCluster*, PeptideHit*> IDPointer;
    typedef boost::variant<const ProteinHit*, const ProteinGroup*, const PeptideCluster*, const PeptideHit*> IDPointerConst;
    //TODO check the impact of different datastructures to store nodes/edges (maybe also use directed graph?)
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer> Graph;
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointerConst> GraphConst;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
    typedef boost::filtered_graph<Graph, boost::function<bool(edge_t)>, boost::function<bool(vertex_t)> > FilteredGraph;


    /// Constructor
    IDBoostGraph(ProteinIdentification &proteins, std::vector<PeptideIdentification>& idedSpectra);

    /// Do sth on connected components
    void applyFunctorOnCCs(std::function<void(FilteredGraph &)> functor);
    void annotateIndistinguishableGroups();


    /// Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type creates a label
    class LabelVisitor:
        public boost::static_visitor<OpenMS::String>
    {
    public:

      OpenMS::String operator()(const PeptideHit* pep) const
      {
        return pep->getSequence().toUnmodifiedString();
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

    };

    /// Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type prints the address. For debug
    class PrintAddressVisitor:
        public boost::static_visitor<>
    {
    public:

      void operator()(PeptideHit* pep) const
      {
        std::cout << pep->getSequence().toUnmodifiedString() << ": " << pep << std::endl;
      }

      void operator()(ProteinHit* prot) const
      {
        std::cout << prot->getAccession() << ": " << prot << std::endl;
      }

      void operator()(const ProteinGroup* protgrp) const
      {
        std::cout << "PG" << std::endl;
      }

      void operator()(const PeptideCluster* /*pc*/) const
      {
        std::cout << "PepClust" << std::endl;
      }

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
        std::cout << "set score " << posterior << " for " << prot->getAccession() << std::endl;
        prot->setScore(posterior);
      }

      void operator()(ProteinGroup* /*protgrp*/, double /*posterior*/) const
      {
        // do nothing
        //protgrp->probability = posterior;
      }

      void operator()(const PeptideCluster* /*pc*/, double /*posterior*/) const
      {
        // do nothing
      }

    };

    /// Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type creates a random
    /// variable or "Dependency" to add into the InferenceGraph
    /*class DependencyVisitor:
        public boost::static_visitor<Dependency>
    {
    public:
      const MessagePasserFactory& mpf;

      DependencyVisitor(const MessagePasserFactory& mpf):
      mpf(mpf)
      {}

      Dependency operator()(const PeptideHit* pep, const std::vector<vertex_t>& neighbors, const FilteredGraph& fg) const
      {
        std::vector<IDBoostGraph::vertex_t> incoming{};
        for (const auto& nb : neighbors)
        {
          if (fg[nb].which() <= 2)
          {
            incoming.push_back(nb);
          }
        }
        if (incoming.size() != 1)
        {
          std::cerr << "Incoming nodes for pep are more than 1. Sth went wrong." << std::endl;
        }
        return mpf.createPeptideEvidenceFactor()
      }

      Dependency operator()(const ProteinHit* prot, const std::vector<vertex_t>& neighbors, const FilteredGraph& fg) const
      {
        return prot->getAccession();
      }

      Dependency operator()(const PeptideCluster* pc, const std::vector<vertex_t>& neighbors, const FilteredGraph& fg) const
      {
        return pep->getSequence().toUnmodifiedString();
      }

      Dependency operator()(const ProteinGroup* pg, const std::vector<vertex_t>& neighbors, const FilteredGraph& fg) const
      {
        return prot->getAccession();
      }

    };*/

    /// Compute connected component on the static graph. Needs to be recomputed if graph is changed.
    void computeConnectedComponents();

    /// Initialize and store the graph
    /// IMPORTANT: Once the graph is built, editing members like (protein/peptide)_hits_ will invalidate it!
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param idedSpectra vector of ProteinIdentifications with links to the proteins and PSMs in its PeptideHits
    /// @param use_all_psms If all or just the FIRST psm should be used
    void buildGraph(bool use_all_psms);
    //void buildGraph(const ProteinIdentification& protein, const std::vector<PeptideIdentification>& peptides);

  private:
    Graph g;
    static PeptideCluster staticPC;
    static ProteinGroup staticPG;
    //GraphConst gconst;
    ProteinIdentification& proteins_;
    std::vector<PeptideIdentification>& idedSpectra_;
    std::vector<unsigned int> componentProperty_;
    unsigned int numCCs_ = 0;

    vertex_t addVertexWithLookup_(IDPointer& ptr, std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map);
    //vertex_t addVertexWithLookup_(IDPointerConst& ptr, std::unordered_map<IDPointerConst, vertex_t, boost::hash<IDPointerConst>>& vertex_map);
  };



} //namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDBOOSTGRAPH_H
