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

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <vector>
#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/variant.hpp>
#include <boost/variant/detail/hash_variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/transform_value_property_map.hpp>

namespace OpenMS
{

  /**
   @brief Creates and maintains a boost graph based on the OpenMS ID datastructures

   For finding connected components.

   @ingroup Analysis_ID
   */
  class IDBoostGraph
  {

  public:

    //typedefs
    typedef boost::variant<const PeptideHit*, const ProteinHit*> IDPointer;
    typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, IDPointer> Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

    /// Constructor
    IDBoostGraph() = default;

    /// Do sth on ccs
    void doSomethingOnCC(const ProteinIdentification& protein, const std::vector<PeptideIdentification>& peptides);


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

    };

    /// Visits nodes in the boost graph (ptrs to an ID Object) and depending on their type creates a random
    /// variable or "Dependency" to add into the InferenceGraph
/*    class RVVisitor:
        public boost::static_visitor<Dependency>
    {
    public:

      Dependency operator()(const PeptideHit* pep, const std::vector<vertex_t>& neighbors) const
      {
        return pep->getSequence().toUnmodifiedString();
      }

      Dependency operator()(const ProteinHit* prot, const std::vector<vertex_t>& neighbors) const
      {
        return prot->getAccession();
      }

    };*/

  private:
    Graph g;
    std::vector<unsigned int> componentProperty;
    unsigned int numCCs = 0;

    void computeConnectedComponents_();

    /// Initialize and store the graph (= maps)
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param peptides vector of ProteinIdentifications with links to the proteins and PSMs in its PeptideHits
    void buildGraph_(const ProteinIdentification& protein, const std::vector<PeptideIdentification>& peptides);

    vertex_t addVertexWithLookup_(const IDPointer& ptr, std::unordered_map<IDPointer, vertex_t, boost::hash<IDPointer>>& vertex_map);

  };



} //namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDBOOSTGRAPH_H
