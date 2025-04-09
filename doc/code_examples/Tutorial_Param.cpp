// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  Param param;

  param.setValue("file:name", "test.xml");
  param.setValue("file:size(MB)", 572.3);
  param.setValue("file:data:min_int", 0);
  param.setValue("file:data:max_int", 16459);

  cout << "Name   : " << (string)(param.getValue("file:name")) << endl;
  cout << "Size   : " << (float)(param.getValue("file:size(MB)")) << endl;
  cout << "Min int: " << (UInt)(param.getValue("file:data:min_int")) << endl;
  cout << "Max int: " << (UInt)(param.getValue("file:data:max_int")) << endl;

  return 0;
} //end of main
