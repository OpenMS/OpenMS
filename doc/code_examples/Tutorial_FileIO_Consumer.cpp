// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/SYSTEM/File.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

class TICWritingConsumer : public MSDataWritingConsumer 
{
  // Inheriting from MSDataWritingConsumer allows to change the data before
  // they are written to disk (to "filename") using the processSpectrum_ and
  // processChromatogram_ functions.
public:
  double TIC {};
  int nr_spectra {};

  // Create new consumer
  TICWritingConsumer(const String& filename) : MSDataWritingConsumer(filename) 
  {}

  // Add a data processing step for spectra before they are written to disk
  void processSpectrum_(MSDataWritingConsumer::SpectrumType & s) override
  {
    for (const auto& p : s) TIC += p.getIntensity();
    nr_spectra++;
  }
  // Empty chromatogram data processing
  void processChromatogram_(MSDataWritingConsumer::ChromatogramType& /* c */) override {}
};

int main(int argc, const char** argv)
{
  // path to the data should be given on the command line
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " <path to tutorial .cpp's, e.g. c:/dev/OpenMS/doc/code_examples/>\n\n";
    return 1;
  }
  String tutorial_data_path(argv[1]);
  auto file_mzXML = tutorial_data_path + "/data/Tutorial_FileIO.mzXML";

  if (! File::exists(file_mzXML)) { std::cerr << "The file " << file_mzXML << " was not found. Did you provide the correct path?\n"; }
  
  // Create the consumer, set output file name, transform
  TICWritingConsumer consumer("Tutorial_FileIO_output.mzML");
  MzMLFile().transform(tutorial_data_path + "/data/Tutorial_FileIO_indexed.mzML", &consumer);

  std::cout << "There are " << consumer.nr_spectra << " spectra in the input file.\n";
  std::cout << "The total ion current is " << consumer.TIC << std::endl;

} //end of main
