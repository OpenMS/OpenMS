// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ScanWindow.h>

using namespace std;

namespace OpenMS
{

  bool ScanWindow::operator==(const ScanWindow & source) const
  {
    return MetaInfoInterface::operator==(source) &&
           begin == source.begin &&
           end == source.end;
  }

  bool ScanWindow::operator!=(const ScanWindow & source) const
  {
    return !(operator==(source));
  }

}

