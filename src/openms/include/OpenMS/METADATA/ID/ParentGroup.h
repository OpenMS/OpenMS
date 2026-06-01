// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/ParentSequence.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief: Group of ambiguously identified parent sequences (e.g. protein group)
    */
    // @TODO: derive from MetaInfoInterface?
    struct ParentGroup
    {
      std::map<ScoreTypeRef, double> scores;
      // @TODO: does this need a "leader" or some such?
      std::set<ParentSequenceRef> parent_refs;
    };

    typedef boost::multi_index_container<
      ParentGroup,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
        boost::multi_index::member<
          ParentGroup, std::set<ParentSequenceRef>,
          &ParentGroup::parent_refs>>>
      > ParentGroups;
    typedef IteratorWrapper<ParentGroups::iterator> ParentGroupRef;

    /** @brief Set of groups of ambiguously identified parent sequences (e.g. results of running a protein inference algorithm)
    */
    struct ParentGroupSet: public ScoredProcessingResult
    {
      String label; // @TODO: use "label" as a uniqueness constraint?
      ParentGroups groups;

      explicit ParentGroupSet(
        const String& label = "",
        const ParentGroups& groups = ParentGroups()):
        label(label), groups(groups)
      {
      }
    };

    typedef std::vector<ParentGroupSet> ParentGroupSets;

  }
}
