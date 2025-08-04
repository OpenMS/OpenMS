// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotatedHit.h>

#include <set>
#include <map>
#include <vector>
#include <algorithm>

namespace OpenMS
{
class OPENMS_DLLAPI NuXLConstants
{
  public:
  static constexpr size_t IA_CHARGE_INDEX = 0;
  static constexpr size_t IA_RANK_INDEX = 1;
  static constexpr size_t IA_DENOVO_TAG_INDEX = 2;
};
}
