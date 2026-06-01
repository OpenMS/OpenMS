// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/METADATA/ID/ParentMatch.h>
#include <OpenMS/METADATA/ID/ScoredProcessingResult.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Representation of an identified sequence (peptide or oligonucleotide)
    template <typename SeqType>
    struct IdentifiedSequence: public ScoredProcessingResult
    {
      SeqType sequence;

      ParentMatches parent_matches;

      explicit IdentifiedSequence(
        const SeqType& sequence,
        const ParentMatches& parent_matches = ParentMatches(),
        const AppliedProcessingSteps& steps_and_scores =
        AppliedProcessingSteps()):
        ScoredProcessingResult(steps_and_scores), sequence(sequence),
        parent_matches(parent_matches)
      {
      }

      IdentifiedSequence(const IdentifiedSequence& other) = default;

      IdentifiedSequence& merge(const IdentifiedSequence& other)
      {
        ScoredProcessingResult::merge(other);
        // merge parent matches:
        for (const auto& pair : other.parent_matches)
        {
          auto pos = parent_matches.find(pair.first);
          if (pos == parent_matches.end()) // new entry
          {
            parent_matches.insert(pair);
          }
          else // merge entries
          {
            pos->second.insert(pair.second.begin(), pair.second.end());
          }
        }

        return *this;
      }

      bool allParentsAreDecoys() const
      {
        if (parent_matches.empty())
        {
          String msg = "no parent found for identified molecule";
          throw Exception::MissingInformation(__FILE__, __LINE__,
                                              OPENMS_PRETTY_FUNCTION, msg);
        }
        for (const auto& pair : parent_matches)
        {
          if (!pair.first->is_decoy) return false;
        }
        return true;
      }
    };

    typedef IdentifiedSequence<AASequence> IdentifiedPeptide;
    typedef IdentifiedSequence<NASequence> IdentifiedOligo;

    // identified peptides indexed by their sequences:
    typedef boost::multi_index_container<
      IdentifiedPeptide,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          IdentifiedPeptide, AASequence, &IdentifiedPeptide::sequence>>>
      > IdentifiedPeptides;
    typedef IteratorWrapper<IdentifiedPeptides::iterator> IdentifiedPeptideRef;

    // identified oligos indexed by their sequences:
    typedef boost::multi_index_container<
      IdentifiedOligo,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          IdentifiedOligo, NASequence, &IdentifiedOligo::sequence>>>
      > IdentifiedOligos;
    typedef IteratorWrapper<IdentifiedOligos::iterator> IdentifiedOligoRef;

  }
}
