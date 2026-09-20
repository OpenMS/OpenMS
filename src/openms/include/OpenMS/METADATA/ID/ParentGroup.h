// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/IDDataContainer.h>

#include <OpenMS/METADATA/ID/ParentSequence.h>


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

    using ParentGroups = IDDataContainer<ParentGroup, std::set<ParentSequenceRef>, std::set<ParentSequenceRef>>;
    extern template class OPENMS_DLLAPI IDDataContainer<ParentGroup, std::set<ParentSequenceRef>, std::set<ParentSequenceRef>>;
    typedef IteratorWrapper<ParentGroups::iterator> ParentGroupRef;

    /** @brief Set of groups of ambiguously identified parent sequences (e.g. results of running a protein inference algorithm)
    */
    struct ParentGroupSet: public ScoredProcessingResult
    {
      std::string label; // @TODO: use "label" as a uniqueness constraint?
      ParentGroups groups;

      explicit ParentGroupSet(
        const std::string& label = "",
        const ParentGroups& groups = ParentGroups()):
        label(label), groups(groups)
      {
      }
    };

    typedef std::vector<ParentGroupSet> ParentGroupSets;

  }
}
