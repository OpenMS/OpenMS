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
      std::string identifier;

      EmpiricalFormula formula;

      std::string name;

      std::string smile;

      std::string inchi;

      explicit IdentifiedCompound(
        const std::string& identifier,
        const EmpiricalFormula& formula = EmpiricalFormula(),
        const std::string& name = "", const std::string& smile = "",
        const std::string& inchi = "", const AppliedProcessingSteps&
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
          IdentifiedCompound, std::string, &IdentifiedCompound::identifier>>>
      > IdentifiedCompounds;
    typedef IteratorWrapper<IdentifiedCompounds::iterator> IdentifiedCompoundRef;
  }
}
