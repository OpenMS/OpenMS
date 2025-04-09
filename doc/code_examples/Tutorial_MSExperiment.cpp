// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_MSExperiment]

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{

  // create a peak map containing 4 dummy spectra and peaks
  MSExperiment exp;

  // The following examples creates a MSExperiment containing four MSSpectrum instances.
  for (Size i = 0; i < 4; ++i)
  {
    MSSpectrum spectrum;
    spectrum.setRT(i);
    spectrum.setMSLevel(1);
    for (float mz = 500.0; mz <= 900; mz += 100.0)
    {
      Peak1D peak;
      peak.setMZ(mz + i);
      spectrum.push_back(peak);
    }

    exp.addSpectrum(spectrum);
  }

  // Iteration over the RT range (2,3) and the m/z range (603,802) and print the peak positions.
  for (auto it = exp.areaBegin(2.0, 3.0, 603.0, 802.0); it != exp.areaEnd(); ++it)
  {
    cout << it.getRT() << " - " << it->getMZ() << endl;
  }

  // Iteration over all peaks in the experiment.
  // Output: RT, m/z, and intensity
  // Note that the retention time is stored in the spectrum (not in the peak object)
  for (auto s_it = exp.begin(); s_it != exp.end(); ++s_it)
  {
    for (auto p_it = s_it->begin(); p_it != s_it->end(); ++p_it)
    {
      cout << s_it->getRT() << " - " << p_it->getMZ() << " " << p_it->getIntensity() << endl;
    }
  }

  // update the data ranges for all dimensions (RT, m/z, int, IM) and print them:
  exp.updateRanges();
  std::cout << "Data ranges:\n";
  exp.printRange(std::cout);
  std::cout << "\nGet maximum intensity on its own: " << exp.getMinMobility() << '\n';
  exp.getMinRT();

  // Store the spectra to a mzML file with:
  FileHandler fh;
  auto tmp_filename = File::getTemporaryFile();
  fh.storeExperiment(tmp_filename, exp, {FileTypes::MZML});

  // And load it with
  fh.loadExperiment(tmp_filename, exp);
  // If we wanted to load only the MS2 spectra we could speed up reading by setting:
  fh.getOptions().setMSLevels({2});
  // and then load from disk: 
  fh.loadExperiment(tmp_filename, exp);

  // note: the file in 'tmp_filename' will be automatically deleted
  return 0;
} // end of main

//! [doxygen_snippet_MSExperiment]
