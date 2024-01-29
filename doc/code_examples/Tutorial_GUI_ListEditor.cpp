// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <QApplication>
#include <OpenMS/VISUAL/ListEditor.h>

using namespace OpenMS;

int main(int argc, char ** argv)
{
  QApplication app(argc, argv);
  ListEditor * listeditor = new ListEditor;
  listeditor->show();
  return app.exec();
}

