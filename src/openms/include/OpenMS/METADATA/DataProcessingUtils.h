// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/DataProcessing.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Queries on data-processing provenance records.

    These queries do not require a feature map or a quality-control metric.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI DataProcessingUtils
  {
  public:
    /**
      @brief Whether a processing record names IsobaricAnalyzer as its software.

      This preserves the historical QCBase/ConsensusMap convention. It does not
      attempt to identify every kind of labeled experiment.

      @param[in] processing Processing records to inspect (an empty list returns false).
    */
    static bool hasIsobaricAnalyzer(const std::vector<DataProcessing>& processing);
  };
} // namespace OpenMS
