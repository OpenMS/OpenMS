// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Init.h>

#include <libxml/parser.h>
#include <mutex>

namespace OpenMS::Internal
{
  void libxmlInitialize()
  {
    static std::once_flag initialized;
    std::call_once(initialized, [] { xmlInitParser(); });
  }
}
