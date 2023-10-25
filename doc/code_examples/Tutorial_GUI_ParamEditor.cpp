// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <QApplication>

using namespace OpenMS;
using namespace std;

Int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  QApplication app(argc, const_cast<char**>(argv));

  Param param;
  ParamXMLFile paramFile;

  paramFile.load(tutorial_data_path + "/data/Tutorial_ParamEditor.ini", param);

  ParamEditor* editor = new ParamEditor(nullptr);
  editor->load(param);
  editor->show();

  app.exec();

  editor->store();
  paramFile.store("Tutorial_ParamEditor_out.ini", param);

  return 0;
} //end of main
