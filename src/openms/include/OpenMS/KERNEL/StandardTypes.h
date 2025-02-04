// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{
  class MSSpectrum;
  class MSChromatogram;
  class MSExperiment;

  //@{
  /**
      @brief Spectrum consisting of raw data points or peaks.

      Meta information includes retention time and MS-level.

      @ingroup Kernel
  */

  typedef MSSpectrum PeakSpectrum;
  /**
      @brief Two-dimensional map of raw data points or peaks.

      @ingroup Kernel
  */
  typedef MSExperiment PeakMap;

  /**
      @brief Chromatogram consisting of raw data points or peaks

      @ingroup Kernel
  */
  typedef MSChromatogram Chromatogram;
  //@}

}


