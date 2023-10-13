// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoDescription.h>

namespace OpenMS
{
  namespace DataArrays
  {

    /// Float data array class
    class FloatDataArray :
      public MetaInfoDescription,
      public std::vector<float>
    {
    };

    /// Integer data array class
    class IntegerDataArray :
      public MetaInfoDescription,
      public std::vector<Int>
    {
    };

    /// String data array class
    class StringDataArray :
      public MetaInfoDescription,
      public std::vector<String>
    {
    };

  }
} // namespace OpenMS

