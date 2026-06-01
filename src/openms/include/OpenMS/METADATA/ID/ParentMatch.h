// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/ParentSequence.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Meta data for the association between an identified molecule (e.g. peptide) and a parent sequence (e.g. protein).
    */
    struct ParentMatch: public MetaInfoInterface
    {
      // in extraordinary cases (e.g. database searches that allow insertions/
      // deletions), the length of the identified molecule may differ from the
      // length of the subsequence in the parent; therefore, store "end_pos":
      Size start_pos, end_pos;

      // String instead of char so modified residues can be represented:
      String left_neighbor, right_neighbor; // neighboring sequence elements

      static constexpr Size UNKNOWN_POSITION = Size(-1);
      static constexpr char UNKNOWN_NEIGHBOR = 'X';
      static constexpr char LEFT_TERMINUS = '[';
      static constexpr char RIGHT_TERMINUS = ']';

      explicit ParentMatch(Size start_pos = UNKNOWN_POSITION,
                                   Size end_pos = UNKNOWN_POSITION,
                                   String left_neighbor = UNKNOWN_NEIGHBOR,
                                   String right_neighbor = UNKNOWN_NEIGHBOR):
        start_pos(start_pos), end_pos(end_pos), left_neighbor(left_neighbor),
        right_neighbor(right_neighbor)
      {
      }

      bool operator<(const ParentMatch& other) const
      {
        // positions determine neighbors - no need to compare those:
        return (std::tie(start_pos, end_pos) <
                std::tie(other.start_pos, other.end_pos));
      }

      bool operator==(const ParentMatch& other) const
      {
        // positions determine neighbors - no need to compare those:
        return (std::tie(start_pos, end_pos) ==
                std::tie(other.start_pos, other.end_pos));
      }

      bool hasValidPositions(Size molecule_length = 0, Size parent_length = 0) const
      {
        if ((start_pos == UNKNOWN_POSITION) || (end_pos == UNKNOWN_POSITION))
        {
          return false;
        }
        if (end_pos < start_pos) return false;
        if (molecule_length && (end_pos - start_pos + 1 != molecule_length))
        {
          return false;
        }
        if (parent_length && (end_pos >= parent_length)) return false;
        return true;
      }
    };

    /// mapping: parent sequence -> match information
    typedef std::map<ParentSequenceRef,
                     std::set<ParentMatch>> ParentMatches;

  }
}
