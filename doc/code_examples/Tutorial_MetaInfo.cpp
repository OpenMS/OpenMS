// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  MetaInfoInterface info;

  // insert meta data
  info.setMetaValue("color", String("#ff0000"));
  info.setMetaValue("id", 112131415);

  // access id by index
  UInt id_index = info.metaRegistry().getIndex("id");
  cout << "id   : " << (UInt)(info.getMetaValue(id_index)) << endl;
  // access color by name
  cout << "color: " << (String)(info.getMetaValue("color")) << endl;

  return 0;
} // end of main
