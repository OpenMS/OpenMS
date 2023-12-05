// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAdductDefinition.h>

using namespace std;

namespace OpenMS
{
  bool NuXLFragmentAdductDefinition::operator<(const NuXLFragmentAdductDefinition& other) const
  {
    String fa = formula.toString();
    String fb = other.formula.toString();
    return std::tie(mass, fa, name) < std::tie(other.mass, fb, other.name);
  }

  bool NuXLFragmentAdductDefinition::operator==(const NuXLFragmentAdductDefinition& other) const
  { 
    return std::tie(formula, name) == std::tie(other.formula, other.name); 
  } 
}

