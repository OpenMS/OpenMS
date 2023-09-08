// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  MzXMLFile mzxml;
  MzMLFile mzml;

  // temporary data storage
  PeakMap map;

  // convert MzXML to MzML
  mzxml.load(tutorial_data_path + "/data/Tutorial_FileIO.mzXML", map);
  mzml.store("Tutorial_FileIO.mzML", map);

  return 0;
} //end of main
