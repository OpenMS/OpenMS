// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! [Precursor]

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  
  if (argc < 2) return 1;
  
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  MSExperiment spectra;

  // load mzML from code examples folder
  FileHandler().loadExperiment(tutorial_data_path + "/data/Tutorial_GaussFilter.mzML", spectra);

  // iterate over map and output MS2 precursor information
  for (auto s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
  {
    // we are only interested in MS2 spectra so we skip all other levels
    if (s_it->getMSLevel() != 2) continue;

    // get a reference to the precursor information
    const MSSpectrum& spectrum = *s_it;
    const vector<Precursor>& precursors = spectrum.getPrecursors();

    // size check & throw exception if needed 
    if (precursors.empty()) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, precursors.size());

    // get m/z and intensity of precursor
    double precursor_mz = precursors[0].getMZ();
    float precursor_int = precursors[0].getIntensity();
  
    // retrieve the precursor spectrum (the most recent MS1 spectrum)
    PeakMap::ConstIterator precursor_spectrum = spectra.getPrecursorSpectrum(s_it);
    double precursor_rt = precursor_spectrum->getRT();
  
    // output precursor information
    std::cout << " precursor m/z: " << precursor_mz
              << " intensity: " << precursor_int
              << " retention time (sec.): " << precursor_rt 
              << std::endl;
   }
                                                            
  return 0;
} // end of main

//! [Precursor]
