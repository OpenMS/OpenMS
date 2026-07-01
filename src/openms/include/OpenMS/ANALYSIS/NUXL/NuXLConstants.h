// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
/**
  @brief Index positions of the per-peak @c IntegerDataArrays used by the OpenNuXL pipeline.

  OpenNuXL annotates each @c MSSpectrum with a fixed layout of parallel
  @c IntegerDataArrays: the array at slot @c i carries the data described by
  @c IA_*_INDEX @c == @c i. Code that accesses these slots indexes them through
  these constants so the layout stays consistent between the spectrum builders in
  @c TOPP_OpenNuXL and the alignment / scoring helpers in
  @ref NuXLAnnotateAndLocate.

  @ingroup Analysis_ID
*/
class OPENMS_DLLAPI NuXLConstants
{
  public:
  /// IntegerDataArray slot holding the per-peak charge.
  static constexpr size_t IA_CHARGE_INDEX = 0;
  /// IntegerDataArray slot holding the per-peak intensity rank (named @c "intensity_rank"; @c 0 marks the most intense peak, ascending).
  static constexpr size_t IA_RANK_INDEX = 1;
  /// IntegerDataArray slot holding the length of the longest amino-acid de novo tag found in the spectrum (named @c "longest_tag"; a single-element array broadcast over the whole spectrum). Always @c 0 unless OpenNuXL is built with the @c CALCULATE_LONGEST_TAG macro defined.
  static constexpr size_t IA_DENOVO_TAG_INDEX = 2;
};
}
