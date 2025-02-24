// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Petra Gutenbrunner $
// $Authors: Petra Gutenbrunner $
// --------------------------------------------------------------------------

//! [doxygen_snippet_MSSpectrum]

#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

int main()
{
  // Create spectrum
  MSSpectrum spectrum;
  Peak1D peak;
  for (float mz = 1500.0; mz >= 500; mz -= 100.0)
  {
    peak.setMZ(mz);
    spectrum.push_back(peak);
  }

  // Sort the peaks according to ascending mass-to-charge ratio
  spectrum.sortByPosition();

  // Iterate over spectrum of those peaks between 800 and 1000 Thomson
  for (auto it = spectrum.MZBegin(800.0); it != spectrum.MZEnd(1000.0); ++it)
  {
    cout << it->getMZ() << endl;
  }
 
  // Access a peak by index
  cout << spectrum[1].getMZ() << " " << spectrum[1].getIntensity() << endl;

  // ... and many more
  return 0;
}

//! [doxygen_snippet_MSSpectrum]
