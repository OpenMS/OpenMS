// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{
  /**
    @brief Metafunction to choose among Peak1D respectively Peak2D through a template argument.

    The result is accessible via typedef Type:
      - @c DPeak<1>::Type is @c Peak1D
      - @c DPeak<2>::Type is @c Peak2D

  */
  template <UInt dimensions>
  struct DPeak
  {};

  // We do not want these classes to show up in the docu
  /// @cond HIDDENSTUFF

  template <>
  struct DPeak<1>
  {
    typedef Peak1D Type;
  };

  template <>
  struct DPeak<2>
  {
    typedef Peak2D Type;
  };

  /// @endcond

} // namespace OpenMS

