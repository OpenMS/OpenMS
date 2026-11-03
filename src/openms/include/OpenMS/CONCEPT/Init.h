// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Exception.h>

/**
    @brief Initialization procedures

    Part of the library needs to be initialized to be used properly, these
    procedures need to be run at startup (and preferably only once).

*/

namespace OpenMS::Internal
{
  /**
    @brief Initialize libxml2 once before parsing from concurrent threads.

    OpenMS deliberately never calls xmlCleanupParser: other threads or dependent
    libraries may still own parser contexts when an OpenMS object is destroyed.
  */
  OPENMS_DLLAPI void libxmlInitialize();
}
