// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinModificationSummary.h>

using namespace OpenMS;

bool ProteinModificationSummary::operator==(const ProteinModificationSummary& rhs) const
{
  return AALevelSummary == rhs.AALevelSummary;
}

bool ProteinModificationSummary::Statistics::operator==(const Statistics& rhs) const
{
  return std::tie(count, frequency, FLR, probability) == std::tie(rhs.count, rhs.frequency, rhs.FLR, rhs.probability);
}

