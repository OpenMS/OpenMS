// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

class TICWritingConsumer : public MSDataWritingConsumer 
{
  // Inheriting from MSDataWritingConsumer allows to change the data before
  // they are written to disk (to "filename") using the processSpectrum_ and
  // processChromatogram_ functions.
public:
  double TIC;
  int nr_spectra;

  // Create new consumer, set TIC to zero
  TICWritingConsumer(String filename) : MSDataWritingConsumer(filename) 
    { TIC = 0.0; nr_spectra = 0;}

  // Add a data processing step for spectra before they are written to disk
  void processSpectrum_(MSDataWritingConsumer::SpectrumType & s) override
  {
    for (Size i = 0; i < s.size(); i++) { TIC += s[i].getIntensity(); }
    nr_spectra++;
  }
  // Empty chromatogram data processing
  void processChromatogram_(MSDataWritingConsumer::ChromatogramType& /* c */) override {}
};

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  // Create the consumer, set output file name, transform
  TICWritingConsumer * consumer = new TICWritingConsumer("Tutorial_FileIO_output.mzML");
  MzMLFile().transform(tutorial_data_path + "/data/Tutorial_FileIO_indexed.mzML", consumer);

  std::cout << "There are " << consumer->nr_spectra << " spectra in the input file." << std::endl;
  std::cout << "The total ion current is " << consumer->TIC << std::endl;
  delete consumer;

  return 0;
} //end of main
