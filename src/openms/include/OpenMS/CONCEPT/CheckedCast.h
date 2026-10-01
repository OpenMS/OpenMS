// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <type_traits>
#include <utility>

namespace OpenMS
{
/// Convert an integer without silently wrapping or truncating an out-of-range value.
template<typename Target, typename Source>
Target checkedCast(Source value)
{
  static_assert(std::is_integral_v<Target> && std::is_integral_v<Source>);
  if (! std::in_range<Target>(value)) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
  return static_cast<Target>(value);
}
} // namespace OpenMS
