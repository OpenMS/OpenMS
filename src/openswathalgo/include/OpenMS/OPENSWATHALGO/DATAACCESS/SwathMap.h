// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

namespace OpenSwath
{
  /**
   * @brief Data structure to hold one SWATH map with information about upper /
   * lower isolation window and whether the map is MS1 or MS2.
   */
  struct SwathMap
  {
    OpenSwath::SpectrumAccessPtr sptr;
    double lower;
    double upper;
    double center;
    double imLower;
    double imUpper;
    bool ms1;

    SwathMap() :
      lower(0.0),
      upper(0.0),
      center(0.0),
      imLower(-1),
      imUpper(-1),
      ms1(false)
    {}

    SwathMap(double mz_start, double mz_end, double mz_center, bool is_ms1)
      : lower(mz_start),
        upper(mz_end),
        center(mz_center),
        ms1(is_ms1)
    {}


    SwathMap(double mz_start, double mz_end, double mz_center, double imLower, double imUpper, bool is_ms1)
      : lower(mz_start),
        upper(mz_end),
        center(mz_center),
	imLower(imLower),
	imUpper(imUpper),
        ms1(is_ms1)
    {}


  };

} //end Namespace OpenSwath

