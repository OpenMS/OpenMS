// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Types.h>
#include <clocale>

namespace OpenMS::Internal
{
  const char * OpenMS_locale = setlocale(LC_ALL, "C");  
} // OpenMS // Internal
