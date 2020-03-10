// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/METADATA/ID/MoleculeParentMatch.h>
#include <OpenMS/METADATA/ID/ScoredProcessingResult.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Identified sequence (peptide or oligonucleotide)
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

      IdentifiedSequence& operator+=(const IdentifiedSequence& other)
      {
        ScoredProcessingResult::operator+=(other);
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
