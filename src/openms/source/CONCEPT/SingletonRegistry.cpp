// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/SingletonRegistry.h>

namespace OpenMS
{

  // unique instance of the SingletonRegistry!
  SingletonRegistry * SingletonRegistry::singletonRegistryInstance_ = nullptr;

} // namespace openms
