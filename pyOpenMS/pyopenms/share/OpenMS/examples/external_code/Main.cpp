// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include "ExampleLibraryFile.h"

using namespace OpenMS;
using namespace OpenMSExternal;

int main(int argc, char * argv[])
{

  FeatureMap fm;
  Feature feature;
  fm.push_back(feature);
  std::string s = ExampleLibraryFile::printSomething();
  std::cout << "From external lib: " << s << "\n";
  std::cout << "All good and well!\n";

  return 0;
}
