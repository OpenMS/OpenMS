// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{
  /**
    @brief Base class for Factory<T>

        Just be able to use dynamic_cast on a pointer

        @ingroup Concept
    */
  class OPENMS_DLLAPI FactoryBase
  {
public:
    /// destructor
    virtual ~FactoryBase(){}

  };

}
