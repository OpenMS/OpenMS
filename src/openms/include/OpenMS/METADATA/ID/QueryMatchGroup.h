// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/MoleculeQueryMatch.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief: Group of related (co-identified) molecule-query matches

      E.g. for cross-linking data or multiplexed spectra.
    */
    struct QueryMatchGroup: public ScoredProcessingResult
    {
      std::set<QueryMatchRef> query_match_refs;

      bool allSameMolecule() const
      {
        // @TODO: return true or false for the empty set?
        if (query_match_refs.size() <= 1) return true;
        const IdentifiedMoleculeRef ref =
          (*query_match_refs.begin())->identified_molecule_ref;
        for (auto it = ++query_match_refs.begin(); it != query_match_refs.end();
             ++it)
        {
          if (!((*it)->identified_molecule_ref == ref)) return false;
        }
        return true;
      }

      bool allSameQuery() const
      {
        // @TODO: return true or false for the empty set?
        if (query_match_refs.size() <= 1) return true;
        DataQueryRef ref = (*query_match_refs.begin())->data_query_ref;
        for (auto it = ++query_match_refs.begin(); it != query_match_refs.end();
             ++it)
        {
          if ((*it)->data_query_ref != ref) return false;
        }
        return true;
      }

      bool operator==(const QueryMatchGroup rhs) const
      {
        return ((rhs.query_match_refs == query_match_refs) &&
                (rhs.steps_and_scores == steps_and_scores));
      }

      bool operator!=(const QueryMatchGroup& rhs) const
      {
        return !operator==(rhs);
      }
    };

    typedef boost::multi_index_container<
      QueryMatchGroup,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
          boost::multi_index::member<QueryMatchGroup, std::set<QueryMatchRef>,
                                     &QueryMatchGroup::query_match_refs>>>
      > QueryMatchGroups;
    typedef IteratorWrapper<QueryMatchGroups::iterator> MatchGroupRef;

  }
}
