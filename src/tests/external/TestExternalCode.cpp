// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>

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
  FileHandler f;
  String tmpfilename = "tmpfile.mzML";

  f.storeExperiment(tmpfilename,exp, {FileTypes::MZML});
  f.loadExperiment(tmpfilename,exp, {FileTypes::MZML});

  std::cout << "Loading and storing of mzML worked!\n";

  std::cout << "All good and well!\n";
  return 0;
}
