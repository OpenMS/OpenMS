// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

// This translation unit holds the out-of-line definition of
// MSExperiment::calculateTIC(). It is intentionally located in the PROCESSING
// layer (not next to the rest of MSExperiment in KERNEL) because the optional
// RT-resampling step uses LinearResamplerAlign, a PROCESSING-layer template.
// Defining the method here keeps the KERNEL -> PROCESSING dependency out of
// the KERNEL layer (PROCESSING -> KERNEL is a legal forward edge).

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>

namespace OpenMS
{
  /// returns the total ion chromatogram (TIC)
  const MSChromatogram MSExperiment::calculateTIC(float rt_bin_size, UInt ms_level) const
  {
    // The TIC is (re)calculated from the MS spectra with set ms_level (default 1).
    // Even if MSExperiment does not contain a TIC chromatogram explicitly, it can be reported.
    MSChromatogram TIC;
    for (const auto& spec: spectra_)
    {
      if ((spec.getMSLevel() == ms_level) || (ms_level == 0))
      {
        // fill chromatogram
        ChromatogramPeakType peak;
        peak.setRT(spec.getRT());
        peak.setIntensity(spec.calculateTIC());
        TIC.push_back(peak);
      }
    }
    if (rt_bin_size > 0)
    {
      LinearResamplerAlign lra;
      Param param = lra.getParameters();
      param.setValue("spacing", rt_bin_size);
      lra.setParameters(param);
      lra.raster(TIC);
    }
    return TIC;
  }

} // namespace OpenMS
