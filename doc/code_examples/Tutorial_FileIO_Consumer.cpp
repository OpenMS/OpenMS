// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data
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
  auto file_mzXML = OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_FileIO_indexed.mzML");
  
  // Create the consumer, set output file name, transform
  TICWritingConsumer consumer("Tutorial_FileIO_output.mzML");
  MzMLFile().transform(file_mzXML, &consumer);

  std::cout << "There are " << consumer.nr_spectra << " spectra in the input file.\n";
  std::cout << "The total ion current is " << consumer.TIC << std::endl;

} //end of main
