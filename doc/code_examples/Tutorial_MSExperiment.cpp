// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [MSExperiment]

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSExperiment.h>
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

  // We could store the spectra to a mzML file with:
  // FileHandler mzml;
  // mzml.storeExperiment(filename, exp);
  
  // And load it with
  // mzml.loadExperiment(filename, exp);
  // If we wanted to load only the MS2 spectra we could speed up reading by setting:
  // mzml.getOptions().addMSLevel(2);
  // before executing: mzml.loadExperiment(filename, exp);

  return 0;
} //end of main

//! [MSExperiment]
