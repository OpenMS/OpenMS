// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser  $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/RNaseDB.h>

using namespace std;

namespace OpenMS
{
  RNaseDB::RNaseDB():
    DigestionEnzymeDB<DigestionEnzymeRNA, RNaseDB>("CHEMISTRY/Enzymes_RNA.xml")
  {
  }

} // namespace OpenMS
