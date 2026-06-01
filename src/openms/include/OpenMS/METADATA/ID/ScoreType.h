// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/MetaData.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Information about a score type.
    */
    struct ScoreType: public MetaInfoInterface
    {
      CVTerm cv_term; // @TODO: derive from CVTerm instead?

      bool higher_better;

      ScoreType():
        higher_better(true)
      {
      }

      explicit ScoreType(const CVTerm& cv_term, bool higher_better):
        cv_term(cv_term), higher_better(higher_better)
      {
      }

      explicit ScoreType(const String& name, bool higher_better):
        cv_term(), higher_better(higher_better)
      {
        cv_term.setName(name);
      }

      ScoreType(const ScoreType& other) = default;

      // don't include "higher_better" in the comparison:
      bool operator<(const ScoreType& other) const
      {
        // @TODO: implement/use "CVTerm::operator<"?
        return (std::tie(cv_term.getAccession(), cv_term.getName()) <
                std::tie(other.cv_term.getAccession(),
                         other.cv_term.getName()));
      }

      // don't include "higher_better" in the comparison:
      bool operator==(const ScoreType& other) const
      {
        return cv_term == other.cv_term;
      }

      bool isBetterScore(double first, double second) const
      {
        if (higher_better) return first > second;
        return first < second;
      }
    };

    typedef std::set<ScoreType> ScoreTypes;
    typedef IteratorWrapper<ScoreTypes::iterator> ScoreTypeRef;
  }
}
