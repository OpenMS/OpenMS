// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <string>
#include <vector>

namespace OpenMS
{

  /**
    @brief Vector of signed integers.

    @ingroup Datastructures
  */
  typedef std::vector<Int> IntList;

  /**
   @brief Vector of double precision real types.

   @ingroup Datastructures
   */
  typedef std::vector<double> DoubleList;


  /**
   @brief Vector of String.

   @ingroup Datastructures
   */
  typedef std::vector<std::string> StringList;

} // namespace OpenMS
