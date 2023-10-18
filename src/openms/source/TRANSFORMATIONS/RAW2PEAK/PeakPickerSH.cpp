// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerSH.h>

#include <OpenMS/KERNEL/SpectrumHelper.h>

namespace OpenMS
{
  PeakPickerSH::PeakPickerSH() :
    DefaultParamHandler("PeakPickerSH"),
    ProgressLogger()
  {
    defaultsToParam_();
  }

  PeakPickerSH::~PeakPickerSH()
  {
    // FLO: Do not care at the moment
  }

  void PeakPickerSH::pickExperiment(const PeakMap & input, PeakMap & output)
  {
    // make sure that output is clear
    output.clear(true);

    // copy experimental settings
    static_cast<ExperimentalSettings &>(output) = input;

    // resize output with respect to input
    output.resize(input.size());

    std::cout << "Before loop, input size = " << input.size() << std::endl;
    Size progress = 0;
    for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
    {
      copySpectrumMeta(input[scan_idx], output[scan_idx]);
      output[scan_idx].setType(SpectrumSettings::CENTROID);

      if (input[scan_idx].getMSLevel() != 1)
      {
        // When not considering MS2 data (MS2 fragment mass tracing=0), Lukas leaves out
        // the entire scan (instead of just copying it to the output as seen in
        // another plugin).
        // pick(input[scan_idx], output[scan_idx], 4.0);
      }
      else
      {
        // TODO: Read value 4.0 from parameters # PeakPickerSH.cpp
        pick(input[scan_idx], output[scan_idx], 5.0);
      }
      setProgress(++progress);
    }
    std::cout << "After loop" << std::endl;

    endProgress();
  }

}
