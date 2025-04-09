// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResampler.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto file_dta = OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_SavitzkyGolayFilter.dta");
  
  // A DTA file always has exactly one Spectrum, so we get that
  MSSpectrum spectrum;
  // Load the dta file into the spectrum
  FileHandler().loadSpectrum(file_dta, spectrum);

  LinearResampler lr;
  Param param_lr;
  param_lr.setValue("spacing", 0.01);
  lr.setParameters(param_lr);
  lr.raster(spectrum);

  SavitzkyGolayFilter sg;
  Param param_sg;
  param_sg.setValue("frame_length", 21);
  param_sg.setValue("polynomial_order", 3);
  sg.setParameters(param_sg);
  sg.filter(spectrum);

  return 0;
} //end of main
