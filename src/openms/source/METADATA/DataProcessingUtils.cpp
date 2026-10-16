// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DataProcessingUtils.h>

#include <algorithm>

namespace OpenMS
{
  bool DataProcessingUtils::hasIsobaricAnalyzer(const std::vector<DataProcessing>& processing)
  {
    return std::any_of(processing.begin(), processing.end(), [](const DataProcessing& dp)
    {
      return dp.getSoftware().getName() == "IsobaricAnalyzer";
    });
  }
} // namespace OpenMS
