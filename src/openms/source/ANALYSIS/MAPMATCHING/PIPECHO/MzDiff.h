// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "OpenMS/DATASTRUCTURES/Param.h"

namespace OpenMS::PipEcho
{

/******************************************************************************/
/**
 * Clean conversion from PPM to m/z.
 */
class MzDiff
{
public:
  MzDiff(const Param& params):
      mz(params.getValue("distance_MZ:max_difference")),
      is_ppm(params.getValue("distance_MZ:unit") == "ppm") {};

  /// Return the maximum allowed difference in m/z which may be
  /// based on the given m/z value if the stored m/z difference is
  /// actually in PPM.
  double mz_diff(const double relative)
  {
    if (is_ppm) { return relative * 1e-6 * mz; }
    else { return mz; }
  };

private:
  double mz;
  bool is_ppm;
};

} // namespace OpenMS::PipEcho
