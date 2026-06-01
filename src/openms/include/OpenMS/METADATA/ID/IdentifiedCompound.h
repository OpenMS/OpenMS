// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/ScoredProcessingResult.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    struct IdentifiedCompound: public ScoredProcessingResult
    {
      String identifier;

      EmpiricalFormula formula;

      String name;

      String smile;

      String inchi;

      explicit IdentifiedCompound(
        const String& identifier,
        const EmpiricalFormula& formula = EmpiricalFormula(),
        const String& name = "", const String& smile = "",
        const String& inchi = "", const AppliedProcessingSteps&
        steps_and_scores = AppliedProcessingSteps()):
        ScoredProcessingResult(steps_and_scores), identifier(identifier),
        formula(formula), name(name), smile(smile), inchi(inchi)
      {
      }

      IdentifiedCompound(const IdentifiedCompound& other) = default;
    };

    // identified compounds indexed by their identifiers:
    typedef boost::multi_index_container<
      IdentifiedCompound,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          IdentifiedCompound, String, &IdentifiedCompound::identifier>>>
      > IdentifiedCompounds;
    typedef IteratorWrapper<IdentifiedCompounds::iterator> IdentifiedCompoundRef;
  }
}
