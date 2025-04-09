// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DTAFile.h>

using namespace std;

namespace OpenMS
{

  DTAFile::DTAFile() :
      default_ms_level_(2) // set default to MS2
  {
  }

  DTAFile::~DTAFile() = default;

} // namespace OpenMS
