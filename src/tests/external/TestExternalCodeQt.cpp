// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Downstream link-contract check: an external consumer that links OpenMS *and*
// wants Qt must find/link Qt itself -- OpenMS no longer re-exports Qt6 as a transitive dependency.
// This translation unit deliberately mixes an OpenMS type with a Qt type; it only builds because the
// accompanying CMakeLists.txt does its own find_package(Qt6) + links Qt6::Core. If OpenMS still leaked
// Qt6 PUBLICly this would compile even without that explicit find (which is exactly what we want to rule out).

#include <OpenMS/KERNEL/MSExperiment.h>
#include <QtCore/QString>
#include <iostream>

using namespace OpenMS;

int main()
{
  // an OpenMS type (provided by the OpenMS package)
  MSExperiment exp;
  exp.resize(2);

  // a Qt type (provided because THIS project found Qt6 itself, not transitively via OpenMS)
  QString msg = QString("OpenMS MSExperiment with %1 spectra; Qt linked explicitly.").arg(exp.size());
  std::cout << msg.toStdString() << '\n';

  return (exp.size() == 2) ? 0 : 1;
}
