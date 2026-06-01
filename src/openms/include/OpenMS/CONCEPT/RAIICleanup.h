// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <functional>

namespace OpenMS
{
  /**
    @brief Exception-safe way of executing arbitrary code at the end of a scope.

    Just pass in a (capturing) lambda function, which will be called upon destruction of an instance of this class.
             
  */

  class RAIICleanup
  {
  public:
    /// no default CTor; we need a lambda
    RAIICleanup() = delete;

    /// pass in any lambda you like which does the cleanup at the end
    RAIICleanup(std::function<void()> l)
      : l_(l)
    {}

    ~RAIICleanup()
    {
      l_();
    }

  private:
    std::function<void()> l_; ///< called upon destruction
  };

} // namespace OPENMS

