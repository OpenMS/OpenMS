// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/ScoredProcessingResult.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Representation of a parent sequence that is identified only indirectly (e.g. a protein).
    */
    struct ParentSequence: public ScoredProcessingResult
    {
      String accession;

      enum MoleculeType molecule_type;

      // @TODO: if there are modifications in the sequence, "sequence.size()"
      // etc. will be misleading!
      String sequence;

      String description;

      double coverage; ///< sequence coverage as a fraction between 0 and 1

      bool is_decoy;

      explicit ParentSequence(
        const String& accession,
        MoleculeType molecule_type = MoleculeType::PROTEIN,
        const String& sequence = "", const String& description = "",
        double coverage = 0.0, bool is_decoy = false,
        const AppliedProcessingSteps& steps_and_scores = AppliedProcessingSteps()):
        ScoredProcessingResult(steps_and_scores), accession(accession),
        molecule_type(molecule_type), sequence(sequence),
        description(description), coverage(coverage), is_decoy(is_decoy)
      {
      }

      ParentSequence(const ParentSequence&) = default;

      ParentSequence& merge(const ParentSequence& other)
      {
        ScoredProcessingResult::merge(other);
        if (sequence.empty()) 
        {
          sequence = other.sequence;
        } 
        else if (!other.sequence.empty() && sequence != other.sequence) // differ and none is empty
        {
          throw Exception::InvalidValue(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, 
                                        "Trying to overwrite ParentSequence sequence '" + sequence + "' with conflicting value.", 
                                        other.sequence);
        } 

        if (description.empty())
        {
          description = other.description;
        } 
        else if (!other.description.empty() && description != other.description) // differ and none is empty
        {
          throw Exception::InvalidValue(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, 
                                        "Trying to overwrite ParentSequence description '" + description + "' with conflicting value.", 
                                        other.description);
        } 

        if (!is_decoy) is_decoy = other.is_decoy; // believe it when it's set
        // @TODO: what about coverage? (not reliable if we're merging data)

        return *this;
      }
    };

    // parent sequences indexed by their accessions:
    // @TODO: allow querying/iterating over proteins and RNAs separately
    typedef boost::multi_index_container<
      ParentSequence,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          ParentSequence, String, &ParentSequence::accession>>>
      > ParentSequences;
    typedef IteratorWrapper<ParentSequences::iterator> ParentSequenceRef;

  }
}
