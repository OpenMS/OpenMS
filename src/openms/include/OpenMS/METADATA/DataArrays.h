// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
      using std::vector<float>::vector; // to allow for aggregate initialization of FloatDataArray
    };

    /// Integer data array class
    class IntegerDataArray :
      public MetaInfoDescription,
      public std::vector<Int>
    {
      using std::vector<int>::vector; // to allow for aggregate initialization of IntegerDataArray
    };

    /// String data array class
    class StringDataArray :
      public MetaInfoDescription,
      public std::vector<String>
    {
      using std::vector<String>::vector; // to allow for aggregate initialization of StringDataArray
    };

  }
} // namespace OpenMS

