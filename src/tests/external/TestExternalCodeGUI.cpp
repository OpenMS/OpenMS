// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// GUI component contract check: this target links OpenMS::OpenMS_GUI, which the
// OpenMS package provides once find_package(OpenMS COMPONENTS GUI) has found the
// Qt6 modules the GUI library was built against; the build succeeds only when
// that link interface resolves. The program itself uses the core library, which
// the GUI library carries as a transitive (PUBLIC) dependency. No GUI header or
// Qt symbol is used, because a library-only installation ships no GUI headers.

#include <OpenMS/KERNEL/MSExperiment.h>

int main()
{
  OpenMS::MSExperiment exp;
  exp.resize(3);
  return (exp.size() == 3) ? 0 : 1;
}
