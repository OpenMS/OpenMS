// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! [doxygen_snippet_MSChromatogram]

#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/ChromatogramSettings.h>

using namespace OpenMS;
using namespace std;

int main()
{
  // create a chromatogram
  MSChromatogram chromatogram;

  // fill it with metadata information
  chromatogram.setNativeID("transition_300.9_188.0");
  chromatogram.getProduct().setMZ(188.0);
  chromatogram.getPrecursor().setMZ(300.9);

  // fill chromatogram with peaks
  ChromatogramPeak peak;
  peak.setIntensity(1.0);
  for (float rt = 200.0; rt >= 100; rt -= 1.0)
  {
    peak.setRT(rt);
    chromatogram.push_back(peak);
  }

  return 0;
} // end of main

//! [doxygen_snippet_MSChromatogram]
