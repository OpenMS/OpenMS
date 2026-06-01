// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data

#include <iostream>


using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto file_mzXML = OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_FileIO_indexed.mzML");

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

} // end of main
