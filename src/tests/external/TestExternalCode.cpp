// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include "ExampleLibraryFile.h"

using namespace OpenMS;
using namespace OpenMSExternal;

int main(int argc, char * argv[])
{
  std::cout << "Call OpenMS function from ExampleLibraryFile" << std::endl;
  ExampleLibraryFile().loadAndSaveFeatureXML();

  std::string s = ExampleLibraryFile::printSomething();
  std::cout << "From external lib: " << s << "\n";

  PeakMap exp;
  MzMLFile f;
  String tmpfilename = "tmpfile.mzML";

  f.store(tmpfilename,exp);
  f.load(tmpfilename,exp);

  std::cout << "Loading and storing of mzML worked!\n";

  std::cout << "All good and well!\n";
  return 0;
}
