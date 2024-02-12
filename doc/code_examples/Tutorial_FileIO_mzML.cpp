// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/SYSTEM/File.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  // path to the data should be given on the command line
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " <path to tutorial .cpp's, e.g. c:/dev/OpenMS/doc/code_examples/>\n\n";
    return 1;
  }
  String tutorial_data_path(argv[1]);
  auto file_mzXML = tutorial_data_path + "/data/Tutorial_FileIO_indexed.mzML";

  if (! File::exists(file_mzXML)) { std::cerr << "The file " << file_mzXML << " was not found. Did you provide the correct path?\n"; }
  
  IndexedMzMLFileLoader imzml;

  // load data from an indexed MzML file
  OnDiscPeakMap map;
  imzml.load(file_mzXML, map);

  // get the first spectrum in memory, do some constant (non-changing) data processing
  MSSpectrum s = map.getSpectrum(0);
  std::cout << "There are " << map.getNrSpectra() << " spectra in the input file." << std::endl;
  std::cout << "The first spectrum has " << s.size() << " peaks." << std::endl;

  // store the (unmodified) data in a different file
  imzml.store("Tutorial_FileIO_output.mzML", map);

} //end of main
