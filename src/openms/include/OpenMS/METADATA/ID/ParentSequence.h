// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/METADATA/ID/ScoredProcessingResult.h>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Representation of a parent sequence that is identified only indirectly (e.g. a protein).
     */
    struct ParentSequence : public ScoredProcessingResult {
      String accession;

      enum MoleculeType molecule_type;

      // @TODO: if there are modifications in the sequence, "sequence.size()"
      // etc. will be misleading!
      String sequence;

      String description;

      double coverage; ///< sequence coverage as a fraction between 0 and 1

      bool is_decoy;

      explicit ParentSequence(const String& accession, MoleculeType molecule_type = MoleculeType::PROTEIN, const String& sequence = "", const String& description = "", double coverage = 0.0,
                              bool is_decoy = false, const AppliedProcessingSteps& steps_and_scores = AppliedProcessingSteps()) :
          ScoredProcessingResult(steps_and_scores),
          accession(accession), molecule_type(molecule_type), sequence(sequence), description(description), coverage(coverage), is_decoy(is_decoy)
      {
      }

      ParentSequence(const ParentSequence&) = default;

      ParentSequence(); // Only for use with Pyopenms FIXME

      ParentSequence& merge(const ParentSequence& other)
      {
        ScoredProcessingResult::merge(other);
        if (sequence.empty())
        {
          sequence = other.sequence;
        }
        else if (!other.sequence.empty() && sequence != other.sequence) // differ and none is empty
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Trying to overwrite ParentSequence sequence '" + sequence + "' with conflicting value.", other.sequence);
        }

        if (description.empty())
        {
          description = other.description;
        }
        else if (!other.description.empty() && description != other.description) // differ and none is empty
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Trying to overwrite ParentSequence description '" + description + "' with conflicting value.", other.description);
        }

        if (!is_decoy)
          is_decoy = other.is_decoy; // believe it when it's set
        // @TODO: what about coverage? (not reliable if we're merging data)

        return *this;
      }
    };

    // parent sequences indexed by their accessions:
    // @TODO: allow querying/iterating over proteins and RNAs separately
    typedef boost::multi_index_container<ParentSequence,
                                         boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::member<ParentSequence, String, &ParentSequence::accession>>>>
      PSeqs;

    struct ParentSequences : public PSeqs {
      ParentSequences() : PSeqs()
      {
      }
      ParentSequences(const ParentSequences& other) : PSeqs(other)
      {
      }
      ParentSequences(const PSeqs& other) : PSeqs(other)
      {
      }
    };

    typedef IteratorWrapper<ParentSequences::iterator, ParentSequence> PSeqR;

    struct ParentSequenceRef : public PSeqR {
      ParentSequenceRef() : PSeqR()
      {
      }
      ParentSequenceRef(const ParentSequenceRef& other) : PSeqR(other)
      {
      }
      ParentSequenceRef(const PSeqR& other) : PSeqR(other)
      {
      }
      ParentSequenceRef(
        const boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<
          boost::multi_index::detail::null_augment_policy,
          boost::multi_index::detail::index_node_base<OpenMS::IdentificationDataInternal::ParentSequence, std::allocator<OpenMS::IdentificationDataInternal::ParentSequence>>>>& other) :
          PSeqR(other)
      {
      }
      ParentSequenceRef operator=(const ParentSequenceRef& other)
      {
        return PSeqR::operator=(other);
      }
    };


  } // namespace IdentificationDataInternal
} // namespace OpenMS
