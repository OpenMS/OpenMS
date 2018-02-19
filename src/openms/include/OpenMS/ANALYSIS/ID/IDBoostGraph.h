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
#include <set>
#include <boost/variant/static_visitor.hpp>

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
      /// Constructor
      IDBoostGraph();

      /// Initialize and store the graph (= maps)
      /// @param protein ProteinIdentification object storing IDs and groups
      /// @param peptides vector of ProteinIdentifications with links to the proteins and PSMs in its PeptideHits
      void buildGraph(const ProteinIdentification& protein, const std::vector<PeptideIdentification>& peptides);
  };

  /// Visits nodes in the boost graph and depending on their type adds random variables to the
  /// Bayesian network
  class InferenceGraphVisitor
      : public boost::static_visitor<>
  {
  public:

    void operator()(PeptideHit* pep) const
    {
      std::cout << "Visited pep: " << pep->getSequence() << std::endl;
    }

    void operator()(ProteinHit* prot) const
    {
      std::cout << "Visited prot: " << prot->getSequence() << std::endl;
    }

  };

  //TODO remove, currently unused
  struct IDBoostGraphNode{
    int activeMember;

    //TODO add other node types
    union {
      const ProteinHit* protein = nullptr;
      const PeptideHit* peptide = nullptr;
    } IDObjectPtr;

  };

  template<typename T>
  struct IDBoostGraphNodeHash {
    inline size_t operator()(const T* pointer) const {
      auto addr = reinterpret_cast<uintptr_t>(pointer);
      #if SIZE_MAX < UINTPTR_MAX
      /* size_t is not large enough to hold the pointer’s memory address */
      addr %= SIZE_MAX; /* truncate the address so it is small enough to fit in a size_t */
      #endif
      return addr;
    }
  };
} //namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDBOOSTGRAPH_H
