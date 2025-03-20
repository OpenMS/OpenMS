// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/ID/ScoreType.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Information about software used for data processing.

      If the same processing is applied to multiple ID runs, e.g. if multiple files (fractions, replicates) are searched with the same search engine, store the software information only once.
    */
    struct ProcessingSoftware: public Software
    {
      /*!
        List of score types assigned by this software, ranked by importance.

        The "primary" score should be the first in the list.
      */
      // @TODO: make this a "list" for cheap "push_front"?
      std::vector<ScoreTypeRef> assigned_scores;

      explicit ProcessingSoftware(
        const String& name = "", const String& version = "",
        const std::vector<ScoreTypeRef>& assigned_scores = std::vector<ScoreTypeRef>()):
        Software(name, version), assigned_scores(assigned_scores)
      {
      }
    };

    // ordering is done using "operator<" inherited from "Software":
    typedef std::set<ProcessingSoftware> ProcessingSoftwares;
    typedef IteratorWrapper<ProcessingSoftwares::iterator> ProcessingSoftwareRef;

  }
}
