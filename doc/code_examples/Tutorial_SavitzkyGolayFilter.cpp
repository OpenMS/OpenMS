// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
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
  
  // A DTA file always has exactly one Spectrum, so we get that
  MSSpectrum spectrum;
  // Load the dta file into the spectrum
  FileHandler().loadSpectrum(tutorial_data_path + "/data/Tutorial_SavitzkyGolayFilter.dta", spectrum);

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
