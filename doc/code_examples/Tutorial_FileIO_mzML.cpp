// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  IndexedMzMLFileLoader imzml;

  // load data from an indexed MzML file
  OnDiscPeakMap map;
  imzml.load(tutorial_data_path + "/data/Tutorial_FileIO_indexed.mzML", map);

  // get the first spectrum in memory, do some constant (non-changing) data processing
  MSSpectrum s = map.getSpectrum(0);
  std::cout << "There are " << map.getNrSpectra() << " spectra in the input file." << std::endl;
  std::cout << "The first spectrum has " << s.size() << " peaks." << std::endl;

  // store the (unmodified) data in a different file
  imzml.store("Tutorial_FileIO_output.mzML", map);

  return 0;
} //end of main
