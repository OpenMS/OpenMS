// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  PeakMap exp_raw;
  PeakMap exp_picked;

  FileHandler().loadExperiment(tutorial_data_path + "/data/Tutorial_PeakPickerCWT.mzML", exp_raw);

  PeakPickerCWT pp;
  Param param;
  param.setValue("peak_width", 0.1);
  pp.setParameters(param);

  pp.pickExperiment(exp_raw, exp_picked);
  exp_picked.updateRanges();

  cout << "\nMinimal fwhm of a mass spectrometric peak: " << (double)param.getValue("peak_width")
       << "\n\nNumber of picked peaks " << exp_picked.getSize() << std::endl;

  return 0;
} //end of main
